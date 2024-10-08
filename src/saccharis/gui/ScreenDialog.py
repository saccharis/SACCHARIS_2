# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'screen_dialog.ui'
#
# Created by: PyQt5 UI code generator 5.15.7
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(797, 514)
        self.user_title = QtWidgets.QLabel(Dialog)
        self.user_title.setGeometry(QtCore.QRect(310, 10, 151, 31))
        self.user_title.setAlignment(QtCore.Qt.AlignCenter)
        self.user_title.setObjectName("user_title")
        self.intersect_title = QtWidgets.QLabel(Dialog)
        self.intersect_title.setGeometry(QtCore.QRect(470, 10, 151, 31))
        self.intersect_title.setAlignment(QtCore.Qt.AlignCenter)
        self.intersect_title.setObjectName("intersect_title")
        self.user_listwidget = QtWidgets.QListWidget(Dialog)
        self.user_listwidget.setGeometry(QtCore.QRect(310, 50, 151, 171))
        self.user_listwidget.setDragEnabled(True)
        self.user_listwidget.setDragDropMode(QtWidgets.QAbstractItemView.DragOnly)
        self.user_listwidget.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.user_listwidget.setSelectionRectVisible(True)
        self.user_listwidget.setObjectName("user_listwidget")
        self.save_counts_button = QtWidgets.QPushButton(Dialog)
        self.save_counts_button.setGeometry(QtCore.QRect(310, 225, 151, 25))
        self.save_counts_button.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.save_counts_button.setObjectName("save_counts_button")
        self.save_categories_button = QtWidgets.QPushButton(Dialog)
        self.save_categories_button.setGeometry(QtCore.QRect(630, 225, 151, 25))
        self.save_categories_button.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.save_categories_button.setObjectName("save_categories_button")
        self.intersection_listwidget = QtWidgets.QListWidget(Dialog)
        self.intersection_listwidget.setGeometry(QtCore.QRect(470, 50, 151, 141))
        self.intersection_listwidget.setDragEnabled(True)
        self.intersection_listwidget.setDragDropMode(QtWidgets.QAbstractItemView.DragOnly)
        self.intersection_listwidget.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.intersection_listwidget.setObjectName("intersection_listwidget")
        self.category_listwidget = QtWidgets.QListWidget(Dialog)
        self.category_listwidget.setGeometry(QtCore.QRect(630, 50, 151, 171))
        self.category_listwidget.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.category_listwidget.setObjectName("category_listwidget")
        self.category_title = QtWidgets.QLabel(Dialog)
        self.category_title.setGeometry(QtCore.QRect(630, 10, 151, 31))
        self.category_title.setAlignment(QtCore.Qt.AlignCenter)
        self.category_title.setObjectName("category_title")
        self.queue_listwidget = QtWidgets.QListWidget(Dialog)
        self.queue_listwidget.setGeometry(QtCore.QRect(310, 310, 311, 161))
        self.queue_listwidget.setDragEnabled(True)
        self.queue_listwidget.setDragDropMode(QtWidgets.QAbstractItemView.DragDrop)
        self.queue_listwidget.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.queue_listwidget.setViewMode(QtWidgets.QListView.IconMode)
        self.queue_listwidget.setObjectName("queue_listwidget")
        self.queue_title = QtWidgets.QLabel(Dialog)
        self.queue_title.setGeometry(QtCore.QRect(310, 281, 321, 31))
        self.queue_title.setAlignment(QtCore.Qt.AlignCenter)
        self.queue_title.setObjectName("queue_title")
        self.cancel_pushButton = QtWidgets.QPushButton(Dialog)
        self.cancel_pushButton.setGeometry(QtCore.QRect(720, 480, 61, 23))
        self.cancel_pushButton.setObjectName("cancel_pushButton")
        self.run_pipeline_button = QtWidgets.QPushButton(Dialog)
        self.run_pipeline_button.setGeometry(QtCore.QRect(630, 430, 151, 41))
        self.run_pipeline_button.setObjectName("run_pipeline_button")
        self.save_intersection_button = QtWidgets.QPushButton(Dialog)
        self.save_intersection_button.setGeometry(QtCore.QRect(470, 225, 151, 25))
        self.save_intersection_button.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.save_intersection_button.setFont(font)
        self.save_intersection_button.setObjectName("save_intersection_button")
        self.add_user_button = QtWidgets.QPushButton(Dialog)
        self.add_user_button.setGeometry(QtCore.QRect(310, 255, 151, 25))
        self.add_user_button.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.add_user_button.setObjectName("add_user_button")
        self.add_intersection_button = QtWidgets.QPushButton(Dialog)
        self.add_intersection_button.setGeometry(QtCore.QRect(630, 255, 151, 25))
        self.add_intersection_button.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.add_intersection_button.setObjectName("add_intersection_button")
        self.add_user_intersect_button = QtWidgets.QPushButton(Dialog)
        self.add_user_intersect_button.setGeometry(QtCore.QRect(470, 255, 151, 25))
        self.add_user_intersect_button.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.add_user_intersect_button.setObjectName("add_user_intersect_button")
        self.save_queue_button = QtWidgets.QPushButton(Dialog)
        self.save_queue_button.setGeometry(QtCore.QRect(630, 390, 151, 31))
        self.save_queue_button.setObjectName("save_queue_button")
        self.include_subfamily_checkbox = QtWidgets.QCheckBox(Dialog)
        self.include_subfamily_checkbox.setGeometry(QtCore.QRect(470, 200, 151, 20))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.include_subfamily_checkbox.setFont(font)
        self.include_subfamily_checkbox.setObjectName("include_subfamily_checkbox")
        self.horizontalLayoutWidget = QtWidgets.QWidget(Dialog)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(310, 480, 311, 26))
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.remove_queue_button_2 = QtWidgets.QPushButton(self.horizontalLayoutWidget)
        self.remove_queue_button_2.setObjectName("remove_queue_button_2")
        self.horizontalLayout_2.addWidget(self.remove_queue_button_2)
        self.clear_queue_button = QtWidgets.QPushButton(self.horizontalLayoutWidget)
        self.clear_queue_button.setObjectName("clear_queue_button")
        self.horizontalLayout_2.addWidget(self.clear_queue_button)
        self.file_list_listWidget = QtWidgets.QListWidget(Dialog)
        self.file_list_listWidget.setGeometry(QtCore.QRect(18, 50, 281, 421))
        self.file_list_listWidget.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.file_list_listWidget.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectItems)
        self.file_list_listWidget.setSelectionRectVisible(True)
        self.file_list_listWidget.setObjectName("file_list_listWidget")
        self.file_label = QtWidgets.QLabel(Dialog)
        self.file_label.setGeometry(QtCore.QRect(20, 10, 281, 31))
        self.file_label.setObjectName("file_label")
        self.save_file_summary_button = QtWidgets.QPushButton(Dialog)
        self.save_file_summary_button.setGeometry(QtCore.QRect(50, 480, 221, 25))
        self.save_file_summary_button.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.save_file_summary_button.setObjectName("save_file_summary_button")

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.user_title.setText(_translate("Dialog", "Families in user FASTA"))
        self.intersect_title.setText(_translate("Dialog", "Families in user FASTA\n"
"and selected category"))
        self.save_counts_button.setText(_translate("Dialog", "Export family selection"))
        self.save_categories_button.setText(_translate("Dialog", "Export category counts"))
        self.category_title.setText(_translate("Dialog", "User categories"))
        self.queue_title.setText(_translate("Dialog", "Family queue for pipeline runs"))
        self.cancel_pushButton.setText(_translate("Dialog", "Cancel"))
        self.run_pipeline_button.setText(_translate("Dialog", "Run Pipeline"))
        self.save_intersection_button.setText(_translate("Dialog", "Export intersect selection"))
        self.add_user_button.setText(_translate("Dialog", "Queue user selection"))
        self.add_intersection_button.setText(_translate("Dialog", "Queue active categories"))
        self.add_user_intersect_button.setText(_translate("Dialog", "Queue intersect selection"))
        self.save_queue_button.setText(_translate("Dialog", "Export queue families"))
        self.include_subfamily_checkbox.setText(_translate("Dialog", "Include subfamilies"))
        self.remove_queue_button_2.setText(_translate("Dialog", "Remove from queue"))
        self.clear_queue_button.setText(_translate("Dialog", "Clear queue"))
        self.file_label.setText(_translate("Dialog", "Fasta file(s):"))
        self.save_file_summary_button.setText(_translate("Dialog", "Export summary for all files"))
