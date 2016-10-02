import sys

from PyQt4.QtCore import *
from PyQt4.QtGui import *
import numpy as np
import queue


from Method import runinterface


class MainApplication(QMainWindow):
    """
    This is the main window for Generalized inversion software
    """

    def __init__(self):


        QMainWindow.__init__(self)
        self.setWindowTitle('Home')
        self.setGeometry(50, 50, 500, 500)


        #Main widget
        self.mainwidget = QWidget()


        #Main Layouts
        self.mainlayout2 = QHBoxLayout()
        self.mainlayout3 = QVBoxLayout()
        self.mainlayout = QVBoxLayout()

        #Stacked layout
        # self.stacked_layout = QStackedLayout()
        # self.stacked_layout.addWidget(self.mainwidget)
        # self.stacked_layout.setCurrentIndex(3)

        #Widgets
        #Form line edits
        self.noffledit = QLineEdit()
        self.noffledit.setText('20')
        self.kh = QLineEdit()
        self.kh.setText('0.045')
        self.kv = QLineEdit()
        self.kv.setText('0.025')
        self.fmin = QLineEdit()
        self.fmin.setText('0.4')
        self.fmax = QLineEdit()
        self.fmax.setText('15')
        self.bs = QLineEdit()
        self.bs.setText('3.5')
        self.ps = QLineEdit()
        self.ps.setText('2.8')
        self.rf = QLineEdit()
        self.alphas = QLineEdit()
        self.alphas.setText('5.98')
        self.Inputfolder = QLineEdit()
        self.Outputfolder = QLineEdit()



        #Labels
        self.programheadings = QLabel('Select which part of program to run')
        # self.programheadings.setGeometry(50, 50, 500, 500)

        #Buttons
        self.run = QPushButton('Run')
        self.inputfolder = QPushButton('Input folder')
        self.outputfolder = QPushButton('Output folder')
        self.Exit = QPushButton('Exit')


        #Checkboxes
        self.interpolation = QCheckBox('Interpolation')
        self.interpolation.setChecked(True)
        self.pga = QCheckBox('Peak Ground Acceleration(PGA)')
        self.pga.setChecked(True)
        self.rfamplification = QCheckBox('Refence-site amplification')
        self.rfamplification.setChecked(True)
        self.matrices = QCheckBox('Making matrices')
        # self.matrices.setChecked(True)
        self.Svd = QCheckBox('SVD calculations')
        # self.Svd.setChecked(True)
        self.path = QCheckBox('Path')
        # self.path.setChecked(True)
        self.site = QCheckBox('Site')
        # self.site.setChecked(True)
        self.source = QCheckBox('Source')
        # self.source.setChecked(True)
        self.gridsearch = QCheckBox('Grid-search')
        # self.gridsearch.setChecked(True)
        self.gridsearchplotting = QCheckBox('Grid-search plots')
        # self.gridsearchplotting.setChecked(True)
        self.extraplots = QCheckBox('Extra plots')
        # self.extraplots.setChecked(True)

        # Form layout creation
        self.form = QFormLayout()
        self.form.addRow('Number of frequencies', self.noffledit)
        self.form.addRow('Horizontal Kappa', self.kh)
        self.form.addRow('Vertical Kappa', self.kv)
        self.form.addRow('Minimum frequency', self.fmin)
        self.form.addRow('Maximum frequency', self.fmax)
        self.form.addRow('S-wave velosity', self.bs)
        self.form.addRow('Density near fault', self.ps)
        self.form.addRow('Reference site number/s', self.rf)
        self.form.addRow('P-wave velosity', self.alphas)
        self.form.addRow('Input folder', self.Inputfolder)
        self.form.addRow('Output folder', self.Outputfolder)

        #Adding stuff to main layouts
        self.mainlayout2.addLayout(self.mainlayout)
        self.mainlayout2.addLayout(self.mainlayout3)
        self.mainlayout3.addWidget(self.programheadings)
        self.mainlayout3.addWidget(self.interpolation)
        self.mainlayout3.addWidget(self.pga)
        self.mainlayout3.addWidget(self.rfamplification)
        self.mainlayout3.addWidget(self.matrices)
        self.mainlayout3.addWidget(self.Svd)
        self.mainlayout3.addWidget(self.path)
        self.mainlayout3.addWidget(self.site)
        self.mainlayout3.addWidget(self.source)
        self.mainlayout3.addWidget(self.gridsearch)
        self.mainlayout3.addWidget(self.gridsearchplotting)
        self.mainlayout3.addWidget(self.extraplots)




        self.mainlayout.addLayout(self.form)
        self.mainlayout.addWidget(self.inputfolder)
        self.mainlayout.addWidget(self.outputfolder)
        self.mainlayout.addWidget(self.run)
        self.mainlayout.addWidget(self.Exit)


        self.mainwidget.setLayout(self.mainlayout2)
        self.mainwidget.setWindowTitle('Home')
        self.mainwidget.setGeometry(50, 50, 800, 600)
        self.mainwidget.show()


        #Connections
        self.inputfolder.clicked.connect(self.inputdialog)
        self.outputfolder.clicked.connect(self.outputdialog)
        self.Exit.clicked.connect(self.exitapplication)
        self.run.clicked.connect(self.mainprogramrun)

    def inputdialog(self):

        name = str(QFileDialog.getExistingDirectory(self, 'Select folder for inputs'))
        self.Inputfolder.setText(name + '/')

    def outputdialog(self):

        name = str(QFileDialog.getExistingDirectory(self, 'Select folder for outputs'))
        self.Outputfolder.setText(name + '/')

    def exitapplication(self):

        sys.exit()

    def mainprogramrun(self):

        kh = float(self.kh.text())
        kv = float(self.kv.text())
        fmin = float(self.fmin.text())
        fmax = float(self.fmax.text())
        bs = float(self.bs.text())
        ps = float(self.ps.text())
        rf = self.rf.text()
        rf = rf.split(',')
        if len(rf) == 1:
            rf = int(rf[0])

        n = int(self.noffledit.text())
        alphas = float(self.alphas.text())
        inputpath = str(self.Inputfolder.text())
        outputpath = str(self.Outputfolder.text())
        interpolationstatus = self.interpolation.isChecked()
        pgastatus = self.pga.isChecked()
        rfamplificationstatus = self.rfamplification.isChecked()
        matricesstatus = self.matrices.isChecked()
        svdstatus = self.Svd.isChecked()
        pathstatus = self.path.isChecked()
        sitestatus = self.site.isChecked()
        sourcestatus = self.source.isChecked()
        gristatus = self.gridsearch.isChecked()
        gridplotstatus = self.gridsearchplotting.isChecked()
        extraplotsstatus = self.extraplots.isChecked()

        variables = [kh, kv, fmin, fmax, bs, ps, rf, n, alphas, inputpath, outputpath,
                 interpolationstatus, pgastatus, rfamplificationstatus, matricesstatus,
                 svdstatus, pathstatus, sitestatus, sourcestatus, gristatus, gridplotstatus,
                 extraplotsstatus]
        print variables
        self.generalrun = runinterface.GeneralizedRun()
        self.generalrun.initUi(self, variables)





def main():

    application = QApplication(sys.argv) #create new application
    QCoreApplication.setAttribute(Qt.AA_X11InitThreads)
    application_window = MainApplication()
    # application_window.show()
    # application_window.raise_() #raise instance to top of window stack
    application.exec_() #monitor application for events



if __name__ == '__main__':
    main()
