###############################################################################
#
#
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# License: GPL v3
###############################################################################
import sys


def ask_yes_no(ask_msg: str, yes_msg: str | None, no_msg: str | None):
    if not ask_msg.__contains__("y/n"):
        ask_msg += " (y/n):"
    try:
        answer = input(ask_msg)
    except KeyboardInterrupt:
        sys.exit(2)
    if answer.lower() == "y" or answer.lower() == "yes":
        if yes_msg:
            print('\n' + yes_msg + ' ')
        return True
    elif answer.lower() == "n" or answer.lower() == "no":
        if no_msg:
            print('\n' + no_msg + ' ')
        return False
    else:
        print("Please input either yes/no [y/n]!")
        return ask_yes_no(ask_msg, yes_msg, no_msg)


def ask_continue():
    answer = input("Would you like to continue the current run of SACCHARIS? (y/n):")
    if answer.lower() == "y" or answer.lower() == "yes":
        print("Continuing...")
    elif answer.lower() == "n" or answer.lower() == "no":
        print("Terminating SACCHARIS...")
        sys.exit(3)
    else:
        print("Please input either yes/no [y/n]!")
        ask_continue()
