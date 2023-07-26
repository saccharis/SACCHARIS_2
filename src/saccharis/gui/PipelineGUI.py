###############################################################################
# GUI interface for SACCHARIS
# Core logic is kept out of all the files in the gui subdirectory as much as possible.
# Original author: Alexander Fraser, https://github.com/AlexSCFraser
# License: GPL v3
###############################################################################
import io
import math
import os
import re
import subprocess
import sys
import logging
from types import SimpleNamespace
from io import TextIOWrapper
from typing import IO

import psutil
from PyQt5.QtCore import QEvent, QObject, QThread, pyqtSignal, pyqtSlot
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QApplication, QListWidgetItem
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QFileDialog
from PyQt5.QtWidgets import QCompleter
from PyQt5.QtWidgets import QDialog
from PyQt5.QtWidgets import QInputDialog
from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtWidgets import QLineEdit
from PyQt5.QtWidgets import QWidget
from PyQt5 import Qt, QtGui

from saccharis.gui import UIDesign
from saccharis.gui import CategoryDialog
from saccharis.gui import SettingsDialog
from saccharis.gui import ScreenDialog

from saccharis.Parse_User_Sequences import concatenate_multiple_fasta
from saccharis.ScreenUserFile import extract_families_hmmer
from saccharis.Pipeline import single_pipeline
from saccharis.CLI import get_version
from saccharis.utils.FamilyCategories import get_category_list, load_family_list_from_file
from saccharis.utils.FamilyCategories import get_user_categories
from saccharis.utils.FamilyCategories import get_default_family_categories
from saccharis.utils.FamilyCategories import write_family_files
from saccharis.utils.FamilyCategories import Matcher
from saccharis.utils.FamilyCategories import save_family_iterable_json
from saccharis.utils.AdvancedConfig import get_user_settings, load_from_env, validate_settings, save_to_file, \
    get_default_settings, get_log_folder, get_output_folder
from saccharis.utils.PipelineErrors import UserError, PipelineException, NewUserFile, make_logger

from saccharis.Cazy_Scrape import Mode
from saccharis.Cazy_Scrape import Domain
from saccharis.ChooseAAModel import TreeBuilder

# local trace function which returns itself
# todo: either delete this or make it some kind of debug option, for sure delete it if the WSL segfault is fixed
# def my_tracer(frame, event, arg = None):
#     # extracts frame code
#     code = frame.f_code
#
#     # extracts calling function name
#     func_name = code.co_name
#
#     # extracts the line number
#     line_no = frame.f_lineno
#
#     print(f"A {event} encountered in \
#     {func_name}() at line number {line_no} ")
#
#     return my_tracer
#
# sys.settrace(my_tracer)

logger = make_logger("GUILogger", get_log_folder(), "gui_logs.txt")

sys._excepthook = sys.excepthook  # save original excepthook


# I don't know if this exception hook code is actually fixing the problem with pyqt5 not catching exception hooks.
# I think the bug is because of an incompatible type issue between a tuple and a list and this is supposed to unpack
# the exception data from either, print the full traceback, and call the original excepthook with individual arguments
# as a workaround to the type incompatibility issue.
def exception_hook(exctype, value, traceback):
    print(exctype, value, traceback)  # print exception.
    print("WARNING: PyQt exception hook bug may or may not be caught right now.")
    sys._excepthook(exctype, value, traceback)  # call original excepthook
    sys.exit(1)


sys.excepthook = exception_hook  # overwrite default excepthook


class BadInputException(Exception):
    def __init__(self, msg, detail=None):
        self.msg = msg
        self.detail = detail
        super().__init__()

    def __str__(self):
        return f"{self.msg} - {self.detail}"


class SACCHARISApp(QMainWindow, UIDesign.Ui_MainWindow):
    kill_signal = pyqtSignal()
    wait_signal = pyqtSignal()
    send_text = pyqtSignal(str)
    send_red_text = pyqtSignal(str)

    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)  # Defined in design.py file by QtDesigner, initializes objects defined by QtDesigner file
        self.setWindowTitle(f"SACCHARIS v{get_version()}")
        # set internal variables
        self.thread = None
        self.worker = None
        self.matcher = Matcher()
        self.queue = None
        self.fam_status = {}
        self.fasta_count_dict = {}
        # connect file browser buttons
        self.select_out_folder_button.clicked.connect(self.browse_folder)
        self.out_folder_lineedit.setText(get_output_folder())
        self.add_fasta_button.clicked.connect(self.add_fasta_files)
        self.family_file_pushbutton.clicked.connect(self.browse_fam_file)
        self.remove_input_button.clicked.connect(self.remove_input_item)
        # set thread dropdown
        self.thread_dropdown.clear()
        for num in range(1, os.cpu_count()+1):
            self.thread_dropdown.addItem(str(num))
        self.thread_dropdown.setCurrentIndex(math.ceil(self.thread_dropdown.count()*.75))
        # setup family autocompleter
        completer = QCompleter(get_category_list("all_families"))
        completer.setCaseSensitivity(Qt.Qt.CaseInsensitive)
        self.family_lineedit.setCompleter(completer)
        # setup family category editing
        self.categories = get_user_categories()
        self.family_categories_dropdown.clear()
        for category in self.categories.keys():
            self.family_categories_dropdown.addItem(category)
        self.add_family_category_button.clicked.connect(self.edit_categories)
        # setup settings editing
        self.settings = get_user_settings()
        self.ncbi_key, self.ncbi_email, self.ncbi_tool = load_from_env(gui_object=self, ask_method=ask_user_yes_no,
                                                                       get_method=get_user_str,
                                                                       show_user_method=tell_user)
        # connect buttons to dialog and run triggers
        self.adv_button.clicked.connect(self.edit_settings)
        self.screen_cazome_button.clicked.connect(self.screen_cazome)
        self.run_button.clicked.connect(self.toggle_logic_thread)
        self.run_button.setCheckable(True)

        # connect text redirection streams
        self.send_text.connect(self.update_text_browser)
        self.console_redirect = TextSignalRedirector(self.send_text)
        sys.stdout = self.console_redirect
        self.send_red_text.connect(self.update_text_browser_red)
        self.err_redirect = TextSignalRedirector(self.send_text)
        # self.err_redirect = TextSignalWrapper(sys.stderr.buffer, self.send_red_text)

        sys.stderr = self.err_redirect
        logger.handlers[0].stream = self.err_redirect

    def closeEvent(self, event):
        print("INFO: GUI close event")
        reply = ask_user_yes_no("Are you sure to quit?", None, None, self)

        if reply == True:
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
            logger.handlers[0].stream = sys.__stderr__
            event.accept()
        else:
            event.ignore()

    def browse_folder(self):
        # self.out_folder_label.clear()
        directory = QFileDialog.getExistingDirectory(self, "Pick a folder for SACCHARIS output", self.out_folder_lineedit.text())
        # directory =
        if directory:
            if os.path.isdir(directory):
                self.out_folder_lineedit.setText(directory)
                # print(self.out_folder_lineedit.text())

    def browse_fam_file(self):
        filepath = QFileDialog.getOpenFileName(caption="Select family file", filter="JSON (*.json)",
                                               directory=self.out_folder_lineedit.text())[0]
        if filepath:
            self.family_file_lineedit.setText(filepath)

    def add_fasta_files(self):
        files = QFileDialog.getOpenFileNames(caption="Select FASTA files", filter="FASTA (*.fasta *.fa)",
                                             directory=self.out_folder_lineedit.text())[0]
        if files:
            for file in files:
                if os.path.isfile(file) and len(self.sequence_source_listwidget.findItems(file, Qt.Qt.MatchExactly)) < 1:
                    self.sequence_source_listwidget.addItem(file)

    def remove_input_item(self):
        this_list = self.sequence_source_listwidget
        removed = this_list.takeItem(this_list.row(this_list.currentItem()))
        # removed.deleteLater()
        # docs say the removed object is not deleted from memory automatically, not sure if the python garbage collector
        # will be able to delete the underlying C++ object, which will cause a memory leak. For some reason, the
        # deleteLater() function normally used for manual object deletion in Qt is not availabe on QListWidgetItem

    def edit_categories(self):
        cat_dialog = CategoryDlg(self, categories=self.categories)
        cat_dialog.exec()
        self.categories = get_user_categories()
        self.family_categories_dropdown.clear()
        self.family_categories_dropdown.addItems(self.categories.keys())
        # print(self.categories)

    def edit_settings(self):
        set_dialog = SettingsDlg(self, settings=self.settings)
        set_dialog.exec()
        self.ncbi_key, self.ncbi_email, self.ncbi_tool = load_from_env(gui_object=self, ask_method=ask_user_yes_no,
                                                                       get_method=get_user_str,
                                                                       show_user_method=tell_user)
        self.settings = get_user_settings()
        # print(self.ncbi_key)
        # print(self.ncbi_email)
        # print(self.ncbi_tool)

    def get_run_options(self):
        if self.family_tabwidget.currentIndex() == 0:
            #  single family
            family = self.family_lineedit.text().strip().upper()
            subfamily = self.subfamily_lineedit.text().strip()
            if subfamily == "":
                subfamily = None
            if family == "":
                raise UserError("Please input a family to run the pipeline on!\n",
                                "If you want to screen a FASTA file for families to run the pipeline on, "
                                "choose the explore tab!")
            elif subfamily and not self.matcher.valid_cazy_family(f"{family}_{subfamily}"):
                raise UserError(f"{family}_{subfamily} is not a valid CAZyme family!")
            elif not self.matcher.valid_cazy_family(family):
                raise UserError(f"{family} is not a valid CAZyme family!")
            args = SimpleNamespace(family_list=[f"{family}{f'_{subfamily}' if subfamily else '' }"],
                                   family_category=None,
                                   explore=False)
        elif self.family_tabwidget.currentIndex() == 1:
            #  family file
            fam_list = load_family_list_from_file(self.family_file_lineedit.text())
            args = SimpleNamespace(family_list=fam_list,
                                   family_category=None,
                                   explore=False)
        elif self.family_tabwidget.currentIndex() == 2:
            # family category
            category = self.family_categories_dropdown.currentText()
            if category not in self.categories:
                raise UserError(f"'{category}' not found in user categories. Please choose an existing category, or "
                                f"add a new one to your user category list.")

            args = SimpleNamespace(family_list=None,
                                   family_category=category,
                                   explore=False)
        else:  # self.family_tabwidget.currentIndex() == 3
            # explore
            if self.queue:  # already got fams to run from user
                args = SimpleNamespace(family_list=self.queue,
                                       family_category=None,
                                       explore=False)
            else:  # have to get fams to run from user
                args = SimpleNamespace(family_list=None,
                                       family_category=None,
                                       explore=True)

        if self.all_radiobutton.isChecked():
            args.__setattr__("cazyme_mode", Mode.ALL_CAZYMES)
        elif self.characterized_radiobutton.isChecked():
            args.__setattr__("cazyme_mode", Mode.CHARACTERIZED)
        else:  # self.structure_radiobutton.isChecked()
            args.__setattr__("cazyme_mode", Mode.STRUCTURE)

        domain_val = 0b0
        if self.archaea_checkbox.isChecked():
            domain_val |= Domain.ARCHAEA.value
        if self.bacteria_checkbox.isChecked():
            domain_val |= Domain.BACTERIA.value
        if self.eukaryota_checkbox.isChecked():
            domain_val |= Domain.EUKARYOTA.value
        if self.viruses_checkbox.isChecked():
            domain_val |= Domain.VIRUSES.value
        if self.unclassified_checkbox.isChecked():
            domain_val |= Domain.UNCLASSIFIED.value
        args.__setattr__("domain", domain_val)

        args.__setattr__("output_path", self.out_folder_lineedit.text())
        user_files = []
        for i in range(self.sequence_source_listwidget.count()):
            user_files.append(self.sequence_source_listwidget.item(i).text())
        args.__setattr__("fasta_file_paths", user_files)
        if len(user_files) < 1 or args.explore:
            args.__setattr__("fasta_file", None)
            args.__setattr__("fasta_source_dict", None)
        else:
            user_merged_file, user_merged_dict = concatenate_multiple_fasta(user_files, output_folder=args.output_path)
            args.__setattr__("fasta_file", user_merged_file)
            args.__setattr__("fasta_source_dict", user_merged_dict)

        if self.fasttree_radiobutton.isChecked():
            args.__setattr__("tree_program", TreeBuilder.FASTTREE)
        else:
            args.__setattr__("tree_program", TreeBuilder.RAXML)

        args.__setattr__("rename", self.auto_prepend_headers_checkbox.isChecked())
        args.__setattr__("force_update", self.fresh_run_checkbox.isChecked())
        args.__setattr__("threads", self.thread_dropdown.currentIndex()+1)
        args.__setattr__("get_fragments", self.include_frag_checkbox.isChecked())
        args.__setattr__("prune_seqs", not self.skip_prune_checkbox.isChecked())
        args.__setattr__("settings", self.settings)

        return args

    def screen_cazome(self):
        try:
            args = self.get_run_options()
        except UserError as error:
            tell_user(error.msg)
            return
        if len(args.fasta_file_paths) < 1:
            tell_user("No fasta file to screen. \nAdd a fasta file and then try screening it!")
        else:
            fasta_count_dict = {}
            self.console_output_textBrowser.clear()

            self.set_user_interaction(False)

            self.console_output_textBrowser.clear()
            self.thread = CazomeScreenThread(self.get_run_options(), self.settings)
            # Connect signals and slots
            # self.thread.started.connect(self.worker.run_pipeline)
            # noinspection PyUnresolvedReferences
            self.thread.finished.connect(self.thread.quit)
            # noinspection PyUnresolvedReferences
            self.kill_signal.connect(self.thread.terminate)
            self.kill_signal.connect(lambda: logger.debug("Kill signal sent"))
            # noinspection PyUnresolvedReferences
            # self.worker.finished.connect(self.worker.deleteLater)
            self.thread.finished.connect(self.thread.deleteLater)
            # noinspection PyUnresolvedReferences
            self.thread.finished.connect(self.report_cazome)
            # noinspection PyUnresolvedReferences
            self.thread.send_status_dict.connect(self.update_family_queue)
            # noinspection PyUnresolvedReferences
            self.thread.send_count_dict.connect(self.set_fasta_count_dict)
            # noinspection PyUnresolvedReferences
            self.thread.send_text.connect(self.update_text_browser)
            # noinspection PyUnresolvedReferences
            self.thread.send_dialog.connect(tell_user)
            self.thread.start()

            # disable buttons and connect reenable for completion
            self.set_user_interaction(False)
            self.run_button.setChecked(True)
            self.run_button.setText("Stop Cazome Screen")

            # scr_dialog = ScreenDlg(fasta_count_dict=fasta_count_dict, categories=self.categories,
            #                        out_dir=self.out_folder_lineedit.text())
            # if scr_dialog.exec():
            #     run_queue = scr_dialog.fams_to_run
            #     self.queue = run_queue
            #     print(run_queue)
            #     try:
            #         self.run()
            #     except UserError as error:
            #         if error.detail_msg:
            #             tell_user(error.msg, error.detail_msg)
            #         else:
            #             tell_user(error.msg)

    def set_fasta_count_dict(self, count_dict: dict):
        self.fasta_count_dict = count_dict

    def report_step(self, step):
        if step == 1:
            self.status_tree_label.setStyleSheet("QLabel{\n"
                                                 "color:Grey\n"
                                                 "}")
            self.cazy_status_label.setStyleSheet("QLabel{\n"
                                                 "color:Green\n"
                                                 "}")
        elif step == 3:
            self.cazy_status_label.setStyleSheet("QLabel{\n"
                                                 "color:Grey\n"
                                                 "}")
            self.status_prune_label.setStyleSheet("QLabel{\n"
                                                  "color:Green\n"
                                                  "}")
        elif step == 4:
            self.status_prune_label.setStyleSheet("QLabel{\n"
                                                  "color:Grey\n"
                                                  "}")
            self.status_alignment_label.setStyleSheet("QLabel{\n"
                                                      "color:Green\n"
                                                      "}")
        elif step == 5:
            self.status_alignment_label.setStyleSheet("QLabel{\n"
                                                      "color:Grey\n"
                                                      "}")
            self.status_mutation_label.setStyleSheet("QLabel{\n"
                                                     "color:Green\n"
                                                     "}")
        elif step == 6:
            self.status_mutation_label.setStyleSheet("QLabel{\n"
                                                     "color:Grey\n"
                                                     "}")
            self.status_tree_label.setStyleSheet("QLabel{\n"
                                                 "color:Green\n"
                                                 "}")

    def report_finished(self):

        self.cazy_status_label.setStyleSheet("QLabel{\n"
                                             "color:Grey\n"
                                             "}")
        self.status_prune_label.setStyleSheet("QLabel{\n"
                                              "color:Grey\n"
                                              "}")
        self.status_alignment_label.setStyleSheet("QLabel{\n"
                                                  "color:Grey\n"
                                                  "}")
        self.status_mutation_label.setStyleSheet("QLabel{\n"
                                                 "color:Grey\n"
                                                 "}")
        self.status_tree_label.setStyleSheet("QLabel{\n"
                                             "color:Grey\n"
                                             "}")
        # reenable user interaction
        self.run_button.setText("Run SACCHARIS")
        self.run_button.setChecked(False)
        self.set_user_interaction(True)
        self.queue = None
        tell_user("Analysis complete!")

    def report_cazome(self):
        self.report_finished()
        scr_dialog = ScreenDlg(fasta_count_dict=self.fasta_count_dict, categories=self.categories,
                               out_dir=os.path.join(self.out_folder_lineedit.text(), "cazome_screen"))
        if scr_dialog.exec():
            run_queue = scr_dialog.fams_to_run
            self.queue = run_queue
            try:
                self.run_button.setChecked(True)  # simulates the button having been pressed
                self.toggle_logic_thread()
            except UserError as error:
                if error.detail_msg:
                    tell_user(error.msg, error.detail_msg)
                else:
                    tell_user(error.msg)

    def set_user_interaction(self, boolean: bool):

        # self.run_button.setEnabled(boolean)
        self.add_fasta_button.setEnabled(boolean)
        self.family_tabwidget.setEnabled(boolean)
        self.select_out_folder_button.setEnabled(boolean)
        self.remove_input_button.setEnabled(boolean)
        self.adv_button.setEnabled(boolean)
        self.archaea_checkbox.setEnabled(boolean)
        self.bacteria_checkbox.setEnabled(boolean)
        self.eukaryota_checkbox.setEnabled(boolean)
        self.viruses_checkbox.setEnabled(boolean)
        self.unclassified_checkbox.setEnabled(boolean)
        self.all_radiobutton.setEnabled(boolean)
        self.characterized_radiobutton.setEnabled(boolean)
        self.structure_radiobutton.setEnabled(boolean)
        self.raxml_radiobutton.setEnabled(boolean)
        self.fasttree_radiobutton.setEnabled(boolean)
        self.fresh_run_checkbox.setEnabled(boolean)
        self.skip_prune_checkbox.setEnabled(boolean)
        self.include_frag_checkbox.setEnabled(boolean)
        self.thread_dropdown.setEnabled(boolean)
        self.auto_prepend_headers_checkbox.setEnabled(boolean)
        self.sequence_source_listwidget.setEnabled(boolean)

    def update_family_queue(self, family_status_dict: dict):
        self.remaining_family_listWidget.clear()
        self.fam_status = family_status_dict
        for family in family_status_dict.items():
            item = QListWidgetItem(family[0])
            if family[1] == 0:
                item.setBackground(QColor(0xBC, 0xBC, 0xBC))  # grey, not yet started
            elif family[1] == 1:
                item.setBackground(QColor(0x5A, 0xCF, 0xC9))  # blue, current pipeline run
                item.setText(item.text() + " - In progress")
            elif family[1] == 2:
                item.setBackground(QColor(0x8A, 0xD8, 0x79))  # green, pipeline complete on item
                item.setText(item.text() + " - Done")
            elif family[1] == 3:
                item.setBackground(QColor(0xF3, 0x53, 0x3A))  # red, error during this item
                item.setText(item.text() + " - ERROR")
            self.remaining_family_listWidget.addItem(item)

    def run(self):
        args = self.get_run_options()
        if args.explore:
            self.screen_cazome()
        else:
            self.console_output_textBrowser.clear()

            # # self.thread = QThread()
            # self.worker = PipelineWorker(args)
            # self.worker.moveToThread(self.thread)  # move worker to the new thread
            # # Connect signals and slots
            # self.thread.started.connect(self.worker.run_pipeline)
            # # noinspection PyUnresolvedReferences
            # self.worker.finished.connect(self.thread.quit)
            # # noinspection PyUnresolvedReferences
            # self.kill_signal.connect(self.thread.quit)
            # self.wait_signal.connect(self.thread.wait)
            # # noinspection PyUnresolvedReferences
            # self.worker.finished.connect(self.worker.deleteLater)
            # self.thread.finished.connect(self.thread.deleteLater)
            # # noinspection PyUnresolvedReferences
            # self.worker.progress_step.connect(self.report_step)
            # # noinspection PyUnresolvedReferences
            # self.worker.progress_family.connect(self.update_family_queue)
            # # noinspection PyUnresolvedReferences
            # self.worker.send_text.connect(self.update_text_browser)
            # self.thread.start()

            self.thread = PipelineThread(args)
            # Connect signals and slots
            # self.thread.started.connect(self.worker.run_pipeline)
            # noinspection PyUnresolvedReferences
            self.thread.finished.connect(self.thread.quit)
            # noinspection PyUnresolvedReferences
            self.kill_signal.connect(self.thread.terminate)
            # noinspection PyUnresolvedReferences
            # self.worker.finished.connect(self.worker.deleteLater)
            self.thread.finished.connect(self.thread.deleteLater)
            # noinspection PyUnresolvedReferences
            self.thread.finished.connect(self.report_finished)
            # noinspection PyUnresolvedReferences
            self.thread.progress_step.connect(self.report_step)
            # noinspection PyUnresolvedReferences
            self.thread.progress_family.connect(self.update_family_queue)
            # noinspection PyUnresolvedReferences
            self.thread.send_text.connect(self.update_text_browser)
            # noinspection PyUnresolvedReferences
            self.thread.send_dialog.connect(tell_user)
            self.thread.start()

            # disable buttons and connect reenable for completion
            self.set_user_interaction(False)

    @pyqtSlot(str)
    def update_text_browser_red(self, string):
        self.update_text_browser(string, QtGui.QColor("Red"))

    def update_text_browser(self, string, color: QtGui.QColor = QtGui.QColor("Black")):
        self.console_output_textBrowser.moveCursor(QtGui.QTextCursor.End)
        self.console_output_textBrowser.ensureCursorVisible()
        self.console_output_textBrowser.setTextColor(color)
        self.console_output_textBrowser.insertPlainText(string)

    def toggle_logic_thread(self):
        logger.debug(f"Run button isChecked: {self.run_button.isChecked()}")
        if self.run_button.isChecked():
            # Starting a new run
            try:
                self.run_button.setText("Stop current run")
                self.run()
            except UserError as error:
                self.run_button.setChecked(False)
                self.run_button.setText("Run SACCHARIS")
                if error.detail_msg:
                    tell_user(error.msg, error.detail_msg)
                else:
                    tell_user(error.msg)
        else:
            # terminating a run
            # sys.stdout.detach()
            print("SACCHARIS pipeline terminated early.")
            logger.warning("SACCHARIS pipeline terminated early.")
            sys.stdout = sys.__stdout__
            # sys.stderr = sys.__stderr__

            for child in psutil.Process().children(recursive=False):
                logger.debug(f"process name(before termination): {child.name()}")
                # if child.name() in ["diamond.exe", "diamond", "wsl.exe", "hmmscan"]:
                if re.match("wsl|diamond|hmmscan|muscle|modeltest-ng|fasttree|raxml", child.name()):
                    logger.debug(f"Attempting to terminate {child.name()}")
                    child.terminate()
            if sys.gettrace():
                print("BEFORE logic thread terminate signal")  # only prints during debug
            logger.debug("BEFORE logic thread terminate signal")
            # noinspection PyUnresolvedReferences
            self.kill_signal.emit()
            if sys.gettrace():
                print("AFTER logic thread terminate signal")  # only prints during debug
            logger.debug("AFTER logic thread terminate signal")
            # self.wait_signal.emit()
            # self.thread.terminate()

            for child in psutil.Process().children(recursive=False):
                logger.debug(f"process name(after termination): {child.name()}")
            curr_fam_item = self.remaining_family_listWidget.findItems("In progress", Qt.Qt.MatchContains).pop()
            curr_fam_item.setBackground(QColor(0xF3, 0x53, 0x3A))  # red, error during this item
            curr_fam_item.setText(curr_fam_item.text().split(' ')[0] + " - ERROR")
            self.run_button.setText("Run SACCHARIS")
            self.set_user_interaction(True)


class CategoryDlg(QDialog):
    """Family category editing dialog."""
    def closeEvent(self, event):
        self.load_user_categories()
        event.accept()

    def __init__(self, parent=None, categories=None):
        super().__init__(parent)
        self.categories = categories
        self.matcher = Matcher()
        self.ui = CategoryDialog.Ui_Dialog()
        self.ui.setupUi(self)
        self.setWindowTitle("Edit Family Categories")
        model = Qt.QStandardItemModel()
        max_width = 0
        for index, name in enumerate(self.categories.keys()):
            group = Qt.QStandardItem(name)
            if max_width < len(name) + 2:
                max_width = len(name) + 2
            for family in self.categories[name]:
                group.appendRow(Qt.QStandardItem(family))
            model.appendRow(group)

        self.ui.category_listwidget.addItems(self.categories.keys())
        self.ui.category_listwidget.currentItemChanged.connect(self.set_families)
        self.ui.delete_family_button.clicked.connect(self.delete_family)
        self.ui.delete_category_button.clicked.connect(self.delete_category)
        self.ui.add_family_pushbutton.clicked.connect(self.add_family)
        self.ui.add_category_pushbutton.clicked.connect(self.add_category)
        completer = QCompleter(get_category_list("all_families"))
        completer.setCaseSensitivity(Qt.Qt.CaseInsensitive)
        self.ui.family_name_linedit.setCompleter(completer)
        self.ui.save_buttonbox.button(self.ui.save_buttonbox.RestoreDefaults).clicked.connect(self.restore_defaults)
        self.ui.save_buttonbox.button(self.ui.save_buttonbox.Cancel).clicked.connect(self.load_user_categories)
        self.ui.save_buttonbox.button(self.ui.save_buttonbox.Save).clicked.connect(self.save_categories)

    def set_families(self):
        selected_category = self.ui.category_listwidget.currentItem()
        self.ui.family_listwidget.clear()
        if selected_category:
            fam_list = self.categories[self.ui.category_listwidget.currentItem().text()]
            self.ui.family_listwidget.addItems(fam_list)

    def delete_family(self):
        to_delete = self.ui.family_listwidget.currentItem()
        if to_delete:
            row = self.ui.family_listwidget.row(to_delete)
            removed_item = self.ui.family_listwidget.takeItem(row)
            self.categories[self.ui.category_listwidget.currentItem().text()].remove(removed_item.text())
        else:
            print("No current family to delete")
            # print(self.categories)

    def delete_category(self):
        to_delete = self.ui.category_listwidget.currentItem()
        if to_delete:
            row = self.ui.category_listwidget.row(to_delete)
            removed_item = self.ui.category_listwidget.takeItem(row)
            self.categories.pop(removed_item.text())
            self.set_families()
        else:
            print("No current category to delete")
        # print(self.categories)

    def add_family(self):
        family = self.ui.family_name_linedit.text().upper()
        if self.ui.category_listwidget.currentItem() is None:
            tell_user("No category selected to add family to!")
        elif self.matcher.valid_cazy_family(family):
            if self.ui.family_listwidget.findItems(family, Qt.Qt.MatchExactly):
                tell_user(f"Family {family} is already in the list!")
            else:
                self.ui.family_listwidget.addItem(family)
                self.categories[self.ui.category_listwidget.currentItem().text()].append(family)
        else:
            tell_user(f"{family} is not a valid CAZy family!!")

    def add_category(self):
        cat_name = self.ui.category_name_lineedit.text()
        if self.ui.category_listwidget.findItems(cat_name, Qt.Qt.MatchExactly):
            tell_user(f"Category name {cat_name} is already in the list!")
        else:
            self.ui.category_listwidget.addItem(cat_name)
            self.categories[cat_name] = []

    def restore_defaults(self):
        self.categories = get_default_family_categories()
        self.ui.category_listwidget.clear()
        self.ui.family_listwidget.clear()
        self.ui.category_listwidget.addItems(self.categories.keys())
        # print(self.categories)

    def save_categories(self):
        write_family_files(None, self.categories)

    def load_user_categories(self):
        self.categories = get_user_categories()
        # print("Cancel")


class SettingsDlg(QDialog):
    def closeEvent(self, event):
        event.accept()

    def __init__(self, parent=None, settings=None):
        super().__init__(parent)
        self.settings = settings
        self.ui = SettingsDialog.Ui_AdvancedSettingsDialog()
        self.ui.setupUi(self)
        self.setWindowTitle("Edit settings")
        self.ncbi_key, self.ncbi_email, self.ncbi_tool = load_from_env(gui_object=self, ask_method=ask_user_yes_no,
                                                                       get_method=get_user_str,
                                                                       show_user_method=tell_user)
        self.set_ui_settings(self.settings)
        self.ui.ncbi_api_lineedit.setText(self.ncbi_key)
        self.ui.ncbi_email_lineedit.setText(self.ncbi_email)
        self.ui.save_button_box.button(self.ui.save_button_box.Save).clicked.connect(self.save_settings)
        self.ui.save_button_box.button(self.ui.save_button_box.RestoreDefaults).clicked.connect(self.restore_defaults)

    def set_ui_settings(self, settings_dict):
        self.ui.hmm_cov_lineedit.setText(str(settings_dict["hmm_cov"]))
        self.ui.hmm_evalue_lineEdit.setText(str(settings_dict["hmm_eval"]))
        self.ui.genbank_querysize_lineedit.setText(str(settings_dict["genbank_query_size"]))
        if settings_dict["raxml_command"].__contains__("AVX2"):
            self.ui.avx2_radiobutton.click()
        elif settings_dict["raxml_command"].__contains__("SSE3"):
            self.ui.sse3_radiobutton.click()
        else:
            self.ui.old_proc_radiobutton.click()
        if settings_dict["raxml_command"].__contains__("PTHREADS"):
            self.ui.multithreading_checkbox.setChecked(True)
        else:
            self.ui.multithreading_checkbox.setChecked(False)

    def get_valid_settings_from_ui(self):
        try:
            raxml_command = "raxmlHPC"
            if self.ui.multithreading_checkbox.isChecked():
                raxml_command += "-PTHREADS"
            if self.ui.avx2_radiobutton.isChecked():
                raxml_command += "-AVX2"
            elif self.ui.sse3_radiobutton.isChecked():
                raxml_command += "-SSE3"

            new_settings = {"hmm_cov": float(self.ui.hmm_cov_lineedit.text()),
                            "hmm_eval": float(self.ui.hmm_evalue_lineEdit.text()),
                            "genbank_query_size": int(self.ui.genbank_querysize_lineedit.text()),
                            "raxml_command": raxml_command
                            }
        except ValueError as error:
            raise BadInputException("Invalid setting input", error.args[0]) from error
        new_ncbi_key = self.ui.ncbi_api_lineedit.text()
        new_email = self.ui.ncbi_email_lineedit.text()
        try:
            validate_settings(new_settings)
            return new_settings, new_ncbi_key, new_email
        except UserWarning as error:
            raise BadInputException(error.args[0])

    def restore_defaults(self):
        self.settings = get_default_settings()
        self.set_ui_settings(self.settings)

    def save_settings(self):
        try:
            settings, key, email = self.get_valid_settings_from_ui()
            # print(settings)
            # print(key)
            # print(email)
            validate_settings(settings)
            save_to_file(settings, email, key)
            self.ncbi_key = key
            self.ncbi_email = email
            tell_user("Successfully updated advanced settings!")
        except BadInputException as error:
            tell_user(error.msg, error.detail)
        except UserWarning as error:
            tell_user("ERROR: Invalid settings, did not save changes.", detail_string=f"Details: {error.args[0]}")


class ScreenDlg(QDialog):
    # def closeEvent(self, event):
    #     event.accept()

    # def __init__(self, parent=None, screened_dict=None, categories):
    def __init__(self, fasta_count_dict, categories, out_dir):
        # super().__init__(parent)
        super().__init__()
        self.fasta_count_dict = fasta_count_dict
        self.active_dict = {}
        self.fams_to_run = []
        self.categories = categories
        self.out_dir = out_dir
        self.matcher = Matcher()
        self.filter_obj = FilterFamily()
        self.ui = ScreenDialog.Ui_Dialog()
        self.ui.setupUi(self)
        self.setWindowTitle("Export/Run CAZome")
        for path in fasta_count_dict.keys():
            self.ui.file_list_listWidget.addItem(path)  # .append(f"\"{os.path.basename(path)}\"
        self.ui.file_list_listWidget.selectAll()
        self.ui.file_list_listWidget.itemSelectionChanged.connect(self.update_active_selection)
        self.ui.file_list_listWidget.itemSelectionChanged.connect(self.update_intersect)
        self.update_active_selection()

        self.ui.category_listwidget.addItems(self.categories.keys())

        self.ui.category_listwidget.itemSelectionChanged.connect(self.update_intersect)
        self.ui.include_subfamily_checkbox.clicked.connect(self.update_intersect)

        self.ui.queue_listwidget.viewport().installEventFilter(self.filter_obj)
        # self.ui.queue_listwidget.drop

        self.ui.cancel_pushButton.clicked.connect(self.close)
        self.ui.add_user_button.clicked.connect(self.queue_user_selection)
        self.ui.add_user_intersect_button.clicked.connect(self.queue_intersect_selection)
        self.ui.add_intersection_button.clicked.connect(self.queue_categories)
        self.ui.remove_queue_button_2.clicked.connect(self.remove_queue_selection)
        self.ui.clear_queue_button.clicked.connect(self.ui.queue_listwidget.clear)
        self.ui.save_counts_button.clicked.connect(self.export_user_selection)
        self.ui.save_intersection_button.clicked.connect(self.export_intersect_selection)
        self.ui.save_categories_button.clicked.connect(self.export_selected_categories)
        self.ui.save_queue_button.clicked.connect(self.export_queue)
        self.ui.save_file_summary_button.clicked.connect(self.export_file_summaries)
        self.ui.run_pipeline_button.clicked.connect(self.run_queue)

    def update_active_selection(self):
        self.active_dict = {}

        for fasta_path_item in self.ui.file_list_listWidget.selectedItems():
            for item in self.fasta_count_dict[fasta_path_item.text()].items():
                if item[0] in self.active_dict:
                    self.active_dict[item[0]] += item[1]
                else:
                    self.active_dict[item[0]] = item[1]

        self.ui.user_listwidget.clear()
        self.ui.user_listwidget.addItems([f"{item[0]}: {item[1]} cazymes" for item in self.active_dict.items()])

    def get_intersect(self):
        intersect = []
        for category in self.ui.category_listwidget.selectedItems():
            cat_list = self.categories[category.text()]
            if self.ui.include_subfamily_checkbox.isChecked():
                intersect_cat = [item for item in self.active_dict.items() if item[0].split('_')[0] in cat_list]
            else:
                intersect_cat = [item for item in self.active_dict.items() if item[0] in cat_list]
            intersect += [item for item in intersect_cat if item not in intersect]

        return intersect

    def update_intersect(self):
        self.ui.intersection_listwidget.clear()
        self.ui.intersection_listwidget.addItems([f"{item[0]}: {item[1]} cazymes" for item in self.get_intersect()])

    def queue_categories(self):
        self.queue_list(self.get_intersect())

    def queue_user_selection(self):
        self.queue_list(self.ui.user_listwidget.selectedItems())

    def queue_intersect_selection(self):
        self.queue_list(self.ui.intersection_listwidget.selectedItems())

    def queue_list(self, item_list):
        try:
            for item in item_list:
                if not self.ui.queue_listwidget.findItems(item.text(), Qt.Qt.MatchExactly):
                    self.ui.queue_listwidget.addItem(item.text())
        except AttributeError:
            strings = [f"{item[0]}: {item[1]} cazymes" for item in item_list]
            for item in strings:
                if not self.ui.queue_listwidget.findItems(item, Qt.Qt.MatchExactly):
                    self.ui.queue_listwidget.addItem(item)

    def remove_queue_selection(self):
        for item in self.ui.queue_listwidget.selectedItems():
            self.ui.queue_listwidget.takeItem(self.ui.queue_listwidget.row(item))

    def export_iterable(self, data, ask_file_path):
        file = QFileDialog.getSaveFileName(caption="Save JSON file", filter="JSON (*.json)",
                                           directory=ask_file_path)[0]
        if file == '':
            return
        try:
            save_family_iterable_json(data, file)
        except IOError as error:
            tell_user("Error occurred while saving file", error.args[0])

    def export_user_selection(self):
        found_file = os.path.join(self.out_dir, re.sub(r"\.fa.*", "_families.json",
                                                       os.path.basename(list(self.fasta_count_dict.keys())[0])))
        data = {item.text().split(':')[0]: int(item.text().split(':')[1].strip().split(' ')[0].strip())
                for item in self.ui.user_listwidget.selectedItems()}
        self.export_iterable(data, found_file)

    def export_intersect_selection(self):
        found_file = os.path.join(self.out_dir, re.sub(r"\.fa.*", "_intersect.json",
                                                       os.path.basename(list(self.fasta_count_dict.keys())[0])))
        data = {item.text().split(':')[0]: int(item.text().split(':')[1].strip().split(' ')[0].strip())
                for item in self.ui.intersection_listwidget.selectedItems()}
        self.export_iterable(data, found_file)

    def export_selected_categories(self):
        found_file = os.path.join(self.out_dir, re.sub(r"\.fa.*", "_category_counts.json",
                                  os.path.basename(list(self.fasta_count_dict.keys())[0])))
        categories_to_save = {}

        for cat in self.ui.category_listwidget.selectedItems():
            if cat.text() in self.categories:
                counts = {family: (self.active_dict[family] if family in self.active_dict else 0)
                          for family in self.categories[cat.text()]}
                categories_to_save[cat.text()] = counts

        self.export_iterable(categories_to_save, found_file)

    def export_queue(self):
        found_file = os.path.join(self.out_dir, re.sub(r"\.fa.*", "_queue.json",
                                                       os.path.basename(list(self.fasta_count_dict.keys())[0])))
        data = {item.text().split(':')[0]: int(item.text().split(':')[1].strip().split(' ')[0].strip())
                for item in [self.ui.queue_listwidget.item(x) for x in range(self.ui.queue_listwidget.count())]}
        self.export_iterable(data, found_file)

    def export_file_summaries(self):
        summary_file = os.path.join(self.out_dir, '_'.join(
            ["counts"] + [os.path.splitext(os.path.basename(name))[0] for name in self.fasta_count_dict.keys()] + ["summary.json"])
                                    )
        self.export_iterable(self.fasta_count_dict, summary_file)

    def run_queue(self):
        self.fams_to_run = [item.text().split(':')[0] for item in
                            [self.ui.queue_listwidget.item(x) for x in range(self.ui.queue_listwidget.count())]]
        self.accept()


class FilterFamily(QObject):
    def eventFilter(self, target: 'QWidget', event: 'QEvent') -> bool:

        if event.type() == QEvent.DragEnter:
            event.acceptProposedAction()
            return True
        # print(event.type())
        # print("drop:", QEvent.Type.Drop)
        if event.type() == QEvent.Drop:
            # handle the event
            # ...
            event.ignore()
            for item in event.source().selectedItems():
                # noinspection PyUnresolvedReferences
                if target.parent().findItems(item.text(), Qt.Qt.MatchExactly):
                    # self.ui.queue_listwidget.addItems("text")
                    # return False/
                    return True
                    pass
                # else:
                    # target.parent().addItem(item.text())
                    # modified = True
                    # return True
        return False  # return false for other event types


class CazomeScreenThread(QThread):
    finished = pyqtSignal()
    early_finished = pyqtSignal()
    progress_step = pyqtSignal(int)
    send_count_dict = pyqtSignal(dict)
    send_status_dict = pyqtSignal(dict)
    send_text = pyqtSignal(str)
    # send_red_text = pyqtSignal(str)
    send_dialog = pyqtSignal(str, str)

    def __init__(self, args, settings):
        super().__init__()
        self.args = args
        self.settings = settings
        # self.console_redirect = StdOutRedirector(sys.stdout.buffer, self.send_text)
        self.console_redirect = TextSignalRedirector(self.send_text)
        # self.err_redirect = TextSignalRedirector(self.send_red_text)

    def run(self):
        file_list = self.args.fasta_file_paths
        cazome_folder = os.path.join(self.args.output_path, "cazome_screen")

        file_status = {file: 0 for file in file_list}
        fasta_count_dict = {}
        # sys.stderr = self.err_redirect

        for fasta_file in file_list:
            # if sys.platform.startswith("win") and not sys.gettrace():
            #     test_dict_list = [{"PL9": 12, "PL9_1": 5, "PL9_4": 1, "GH1": 5, "GH5": 7, "GH16": 4, "AA5": 3},
            #                       {"PL9": 12, "PL7_1": 5, "PL7": 1, "GH2": 5, "GH15": 7, "GH26": 4, "AA1": 3},
            #                       {"PL9": 12, "PL8_1": 1, "GH3": 5, "GH25": 7, "GH36": 4, "AA2": 3}]
            #     # noinspection PyUnresolvedReferences
            #     self.send_dialog.emit("Windows support is not yet fully implemented, for GUI functionality use a linux "
            #                           "OS.\nNow loading test data for GUI testing...")
            #     family_dict = test_dict_list.pop()  # todo: remove fake test data from windows
            # else:
            file_status[fasta_file] = 1
            # noinspection PyUnresolvedReferences
            self.send_status_dict.emit(file_status)
            try:
                sys.stdout = self.console_redirect
                family_dict = extract_families_hmmer(fasta_file, cazome_folder, self.args.threads,
                                                     self.settings["hmm_eval"], self.settings["hmm_cov"])
                hmmer_file = os.path.join(cazome_folder, "hmmer.out")
                ren_file = os.path.join(cazome_folder, f"{os.path.splitext(os.path.basename(fasta_file))[0]}_hmmer.out")
                os.replace(hmmer_file, ren_file)
                os.remove(os.path.join(cazome_folder, "overview.txt"))
                os.remove(os.path.join(cazome_folder, "uniInput"))

            except PipelineException as error:
                # noinspection PyUnresolvedReferences
                self.send_dialog.emit(f"Error while analyzing file {fasta_file}", error.msg)
                # noinspection PyUnresolvedReferences
                self.send_text.emit(error.msg)
                file_status[fasta_file] = 3
                family_dict = None
            except UserWarning as error:
                # noinspection PyUnresolvedReferences
                self.send_dialog.emit(error.args[0], error.args[0])
                # noinspection PyUnresolvedReferences
                self.send_text.emit(error.args[0])
                file_status[fasta_file] = 3
                family_dict = None
                break
            except Exception as error:
                logger.exception("Unhandled exception", exc_info=True)
                # noinspection PyUnresolvedReferences
                self.send_dialog.emit("UNHANDLED EXCEPTION (Please report as a bug on github):", str(error))
                # noinspection PyUnresolvedReferences
                self.send_text.emit(f"UNHANDLED EXCEPTION (Please report as a bug on github): {str(error)}")
                file_status[fasta_file] = 3
                family_dict = None
                break
            finally:
                sys.stdout = sys.__stdout__
            # todo: save data into md5 hash files for future speed
            if family_dict:
                fasta_count_dict[fasta_file] = family_dict
                file_status[fasta_file] = 2
            # noinspection PyUnresolvedReferences
            self.send_status_dict.emit(file_status)

        # sys.stderr = sys.__stderr__
        # noinspection PyUnresolvedReferences
        self.send_count_dict.emit(fasta_count_dict)
        # noinspection PyUnresolvedReferences
        self.send_status_dict.emit(file_status)
        # noinspection PyUnresolvedReferences
        self.finished.emit()


class PipelineThread(QThread):
    finished = pyqtSignal()
    early_finished = pyqtSignal()
    progress_step = pyqtSignal(int)
    progress_family = pyqtSignal(dict)
    send_text = pyqtSignal(str)
    # send_red_text = pyqtSignal(str)
    send_dialog = pyqtSignal(str, str)

    def __init__(self, args):
        super().__init__()
        self.args = args
        # self.console_redirect = StdOutRedirector(sys.stdout.buffer, self.send_text)
        self.console_redirect = TextSignalRedirector(self.send_text)
        # self.err_redirect = TextSignalRedirector(self.send_red_text)

    def run(self):
        if self.args.explore:
            # noinspection PyUnresolvedReferences
            self.send_dialog.emit("Something has gone wrong, this code should not execute!", "args.explore erroneously "
                                                                                             "True")
            # noinspection PyUnresolvedReferences
            self.send_text.emit("Something has gone wrong, this code should not execute!")
            fam_list = []
        elif self.args.family_category:
            fam_list = get_category_list(self.args.family_category)
        else:
            fam_list = self.args.family_list

        # sys.stderr = self.err_redirect

        fam_status = {family: 0 for family in fam_list}
        for fam in fam_list:
            fam_status[fam] = 1
            # noinspection PyUnresolvedReferences
            self.progress_family.emit(fam_status)
            sys.stdout = self.console_redirect
            try:
                single_pipeline(fam, self.args.output_path, self.args.cazyme_mode, domain_mode=self.args.domain,
                                threads=self.args.threads, tree_program=self.args.tree_program,
                                get_fragments=self.args.get_fragments, prune_seqs=self.args.prune_seqs, verbose=False,
                                force_update=self.args.force_update, user_file=self.args.fasta_file,
                                auto_rename=self.args.rename, settings=self.args.settings,
                                gui_step_signal=self.progress_step, merged_dict=self.args.fasta_source_dict)
                fam_status[fam] = 2
            except NewUserFile as error:
                self.args.fasta_file = error.msg
                try:
                    single_pipeline(fam, self.args.output_path, self.args.cazyme_mode, domain_mode=self.args.domain,
                                    threads=self.args.threads, tree_program=self.args.tree_program,
                                    get_fragments=self.args.get_fragments, prune_seqs=self.args.prune_seqs,
                                    verbose=False, force_update=self.args.force_update, user_file=self.args.fasta_file,
                                    auto_rename=self.args.rename, settings=self.args.settings,
                                    gui_step_signal=self.progress_step, merged_dict=self.args.fasta_source_dict)
                    fam_status[fam] = 2
                except PipelineException as error:
                    # noinspection PyUnresolvedReferences
                    self.send_dialog.emit(f"Error while analyzing family {fam}", error.msg)
                    # noinspection PyUnresolvedReferences
                    self.send_text.emit(error.msg)
                    fam_status[fam] = 3
                except UserWarning as error:
                    # noinspection PyUnresolvedReferences
                    self.send_dialog.emit(error.args[0], error.args[0])
                    # noinspection PyUnresolvedReferences
                    self.send_text.emit(error.args[0])
                    # noinspection PyUnresolvedReferences
                    break
            except PipelineException as error:
                # noinspection PyUnresolvedReferences
                self.send_dialog.emit(f"Error while analyzing family {fam}", error.msg)
                # noinspection PyUnresolvedReferences
                self.send_text.emit(error.msg)
                fam_status[fam] = 3
            except UserWarning as error:
                # noinspection PyUnresolvedReferences
                self.send_dialog.emit(f"USER WARNING: {error.args[0]}", error.args[0])
                # noinspection PyUnresolvedReferences
                self.send_text.emit(error.args[0])
                fam_status[fam] = 3
                # noinspection PyUnresolvedReferences
                break
            except Exception as error:
                logger.exception("Unhandled exception", exc_info=True)
                # noinspection PyUnresolvedReferences
                self.send_dialog.emit("UNHANDLED EXCEPTION (Please report as a bug on github):", str(error))
                # noinspection PyUnresolvedReferences
                self.send_text.emit(f"UNHANDLED EXCEPTION (Please report as a bug on github): {str(error)}")
                fam_status[fam] = 3
                # return
                break
            finally:
                sys.stdout = sys.__stdout__
            # self.progress_step.emit(i + 1)

        # sys.stderr = sys.__stderr__
        # noinspection PyUnresolvedReferences
        self.progress_family.emit(fam_status)
        # noinspection PyUnresolvedReferences
        self.finished.emit()


class PipelineWorker(QObject):
    # todo: figure out why this class is here. Why am I not using pipelineworker, did i mean to delete this?
    #  did i mean to switch to QObject instead of QThread?
    finished = pyqtSignal()
    progress_step = pyqtSignal(int)
    progress_family = pyqtSignal(dict)
    send_text = pyqtSignal(str)
    send_red_text = pyqtSignal(str)

    def __init__(self, args):
        super().__init__()
        self.args = args
        # self.console_redirect = StdOutRedirector(sys.stdout.buffer, self.send_text)
        self.console_redirect = TextSignalRedirector(self.send_text)
        # self.err_redirect = TextSignalRedirector(self.send_red_text)

    def run_pipeline(self):
        if self.args.explore:
            tell_user("Something has gone wrong, this code should not execute!")
            fam_list = []
        elif self.args.family_category:
            fam_list = get_category_list(self.args.family_category)
        else:
            fam_list = self.args.family_list

        # sys.stderr = self.err_redirect

        fam_status = {family: 0 for family in fam_list}
        for fam in fam_list:
            fam_status[fam] = 1
            # noinspection PyUnresolvedReferences
            self.progress_family.emit(fam_status)
            try:
                sys.stdout = self.console_redirect
                single_pipeline(fam, self.args.output_path, self.args.cazyme_mode, domain_mode=self.args.domain,
                                threads=self.args.threads, tree_program=self.args.tree_program,
                                get_fragments=self.args.get_fragments, prune_seqs=self.args.prune_seqs, verbose=False,
                                force_update=self.args.force_update, user_file=self.args.fasta_file,
                                auto_rename=self.args.rename, settings=self.args.settings,
                                gui_step_signal=self.progress_step, merged_dict=self.args.fasta_source_dict)
                fam_status[fam] = 2
            except NewUserFile as error:
                self.args.fasta_file = error.msg
                try:
                    single_pipeline(fam, self.args.output_path, self.args.cazyme_mode, domain_mode=self.args.domain,
                                    threads=self.args.threads, tree_program=self.args.tree_program,
                                    get_fragments=self.args.get_fragments, prune_seqs=self.args.prune_seqs,
                                    verbose=False, force_update=self.args.force_update, user_file=self.args.fasta_file,
                                    auto_rename=self.args.rename, settings=self.args.settings,
                                    gui_step_signal=self.progress_step, merged_dict=self.args.fasta_source_dict)
                    fam_status[fam] = 2
                except PipelineException as error:
                    fam_status[fam] = 3
                    tell_user(error.msg)
                except UserWarning as error:
                    tell_user(error.args[0])
                    # noinspection PyUnresolvedReferences
                    self.finished.emit()
                    return
            except PipelineException as error:
                fam_status[fam] = 3
                tell_user(error.msg)
            except UserWarning as error:
                tell_user(error.args[0])
                # noinspection PyUnresolvedReferences
                self.finished.emit()
                return
            except Exception as error:
                tell_user("UNHANDLED EXCEPTION (Please report as a bug on github):", error.args[0])
            finally:
                sys.stdout = sys.__stdout__
            # self.progress_step.emit(i + 1)

        # sys.stderr = sys.__stderr__
        # noinspection PyUnresolvedReferences
        self.progress_family.emit(fam_status)
        # noinspection PyUnresolvedReferences
        self.finished.emit()


class TextSignalRedirector(io.StringIO):
    def __init__(self, update_ui: pyqtSignal, subprocess_file_descriptor=None):
        # super().__init__(buffer)
        super().__init__()
        self.update_ui = update_ui
        self.file_descriptor = subprocess_file_descriptor
        self.null_file_descriptor = os.open(os.devnull, os.O_RDWR)

    def write(self, string):
        # noinspection PyUnresolvedReferences
        self.update_ui.emit(string)
        # self.text_browser.moveCursor(QtGui.QTextCursor.End)
        # self.text_browser.ensureCursorVisible()
        # self.text_browser.insertPlainText(string)

        # self.text_browser.append(string)
        # self.console_window.moveCursor(QtGui.QTextCursor.End)
        # self.console_window.ensureCursorVisible()
        # self.console_window.insertPlainText(text)

    def close(self) -> None:
        os.close(self.null_file_descriptor)
        super().close()

    def fileno(self):
        # return a file descriptor passed in when this object was created for subprocesses to write to when the
        #  stream cannot be redirected to gui textbox
        if self.file_descriptor:
            return self.file_descriptor
        else:
            return self.null_file_descriptor


class TextSignalWrapper(io.TextIOWrapper):
    def __init__(self, buffer, update_ui: pyqtSignal):
        # super().__init__(buffer)
        super().__init__(buffer)
        self.update_ui = update_ui

    def write(self, string):
        # noinspection PyUnresolvedReferences
        self.update_ui.emit(string)
        # self.text_browser.moveCursor(QtGui.QTextCursor.End)
        # self.text_browser.ensureCursorVisible()
        # self.text_browser.insertPlainText(string)

        # self.text_browser.append(string)
        # self.console_window.moveCursor(QtGui.QTextCursor.End)
        # self.console_window.ensureCursorVisible()
        # self.console_window.insertPlainText(text)

    def close(self) -> None:
        # this is here so that you can put a breakpoint on the close call for easier debugging
        super().close()


def tell_user(string, detail_string=None):
    msg_box = QMessageBox()
    msg_box.setText(string)
    msg_box.setWindowTitle("Information")
    if detail_string:
        msg_box.setDetailedText(detail_string)
    msg_box.exec()


def get_user_str(title_msg, item_msg, gui_parent):
    text, ok = QInputDialog().getText(gui_parent, title_msg,
                                      item_msg, QLineEdit.Normal)
    if ok and text:
        return text
    else:
        return None


def ask_user_yes_no(question_string, yes_msg, no_msg, gui_parent):
    # msg_box = QMessageBox()
    button_reply = QMessageBox.question(gui_parent, "Answer required", question_string, QMessageBox.Yes, QMessageBox.No)
    # button_reply = QMessageBox.question("Test", question_string, QMessageBox.Yes, QMessageBox.No)

    if button_reply == QMessageBox.Yes:
        if yes_msg:
            QMessageBox.information(gui_parent, "Accepted", yes_msg)
        return True
    else:
        if no_msg:
            QMessageBox.information(gui_parent, "Cancelled", no_msg)
        return False


def main():
    app = QApplication(sys.argv)
    form = SACCHARISApp()
    form.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
