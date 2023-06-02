###############################################################################
#
#
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# License: GPL v3
###############################################################################


class PipelineException(Exception):
    def __init__(self, msg):
        self.msg = msg
        super().__init__()

    def __str__(self):
        return self.msg


class AAModelError(PipelineException):
    def __init__(self, msg):
        super().__init__(msg)


class FileError(PipelineException):
    def __init__(self, msg):
        super().__init__(msg)


class UserError(PipelineException):
    def __init__(self, msg, detail=None):
        super().__init__(msg)
        self.detail_msg = detail


class NewUserFile(PipelineException):
    def __init__(self, msg):
        super().__init__(msg)


