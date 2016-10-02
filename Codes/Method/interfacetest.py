from PyQt4 import QtCore



class A(QtCore.QThread):

    signal = QtCore.pyqtSignal()
    def __init__(self):
        QtCore.QThread.__init__(self)


    def signalemit(self):
        self.signal.emit()

    def call(self):
        self.signalemit()

    def run(self):
        self.call()


class B():

    def __init__(self):
       pass

    def thread(self):
        th = A()
        th.start()
        th.signal.connect(self.printt)

    def printt(self):

        print 'aaaaaa'


g = B()
g.thread()
