import sys

import os
import queue
import time

from PyQt4.QtCore import *
from PyQt4.QtGui import *
import numpy as np

from Method import effectscalculations,\
    GeneralizedInversion,\
    interpolating,\
    rocksiteamplification

class EffectsCalculationThread(QThread):

    signal = pyqtSignal(str, bool)


    def __init__(self, effect, parent):
        QThread.__init__(self)
        self.variables = effect



    def effectmethod(self):

        # print ' effect methods thread is ruunig for part number %d' % self.effect[22]


        effect_counter = 0
        self.variables.append(effect_counter)
        if self.variables[13]:


            referenceamplification, log = rocksiteamplification.BooreGenericRockSite(self.variables[2],
                                                                                self.variables[3])
            self.signalemit('Reference site amplification is running', True)
            self.signalemit(log, False)

        if self.variables[11]:
            self.signalemit('Interpolating data:', True)
            self.freqs, self.spec, self.rat =  interpolating.interpolating(self.variables[9], self.variables[10]
                                                            ,self.variables[7],self.variables[2],
                                                            self.variables[3])
            self.signalemit('Interpolating data is done', False)


        if self.variables[12]:
            self.signalemit('PGA figure is calculating', True)
            effectscalculations.PGA(self.variables[10], self.spec)

            self.signalemit('PGA figure is saved', False)

        ginversion = GeneralizedInversion.MatrixNSVD()
        if self.variables[14]:
            self.variables[22] = 3
            self.signalemit('Constructing Data matrix, Model matrix:', True)
            self.g, self.datamatrix, self.eqcount, self.stcount, log = ginversion.gmatrix(self.freqs, self.spec, self.rat
                                                                           , self.variables[10], self.variables[6],
                                                                           referenceamplification,
                                                                           self.variables[1], self.variables[4],
                                                                           self.variables[7])

            self.signalemit('Matrixes are done constructing.', True)
            self.signalemit(log, False)
        if self.variables[15]:
            self.variables[22] = 4
            self.signalemit('Sigular value decomposition method is executing:', True)
            results, log = ginversion.svdpart(self.variables[10], self.g, self.datamatrix)
            np.savetxt(self.variables[10] + 'svd/answer.txt', results,
                       fmt='%0.5f', delimiter=',')
            self.signalemit(log, True)
            self.signalemit('SVD method is implemented and results-matrix is calculated', False)

        if self.variables[16]:
            self.signalemit('Path effect is calculating using result-matrix:', True)
            effectscalculations.pathcalculations(self.variables[10], self.variables[7])
            self.signalemit('Path effect is done calculating and the plot is saved', False)
        if self.variables[17]:
            self.signalemit('Site effect is calculating:', True)
            self.dk = self.variables[0] - self.variables[1]
            self.HtoV = (False, self.dk)
            effectscalculations.siteextraction(self.variables[9], self.variables[10],
                                               self.eqcount, self.stcount,
                                               self.variables[7], self.HtoV)
            self.signalemit('Site effects are done calculating and plot for each site is saved.', False)
        if self.variables[18]:
            self.signalemit('Source effects are calculating:', True)
            effectscalculations.sourceextraction(self.variables[9], self.variables[10], self.eqcount, self.variables[4]
                                                 , self.variables[5], self.variables[7])
            self.signalemit('Source effects are dont calculating and plots are saved', False)
        if self.variables[19]:
            self.signalemit('Grid-search method is running using predefined range for parameters', True)
            sigma = (1, 600, 600)   # (Start, End, Number of samples)
            gamma = (2.0, 2.0, 1)    # (Start, End, Number of samples)
            magnitudes = (0.5, 40)   # (Magnitude increase limit, Step)
            effectscalculations.gridsearch(self.variables[9], self.variables[10], sigma,
                                           gamma, magnitudes, self.variables[4],
                                           self.variables[7], self.variables[5], self.variables[8])
            self.signalemit('Grid-search method is completed and the results is saved', False)
        if self.variables[20]:
            effectscalculations.gridsearchplots(self.variables[9], self.variables[10], self.eqcount, self.stcount)

        if self.variables[21]:
            self.signalemit('Extra plots for source mechanism is calculating:', True)
            effectscalculations.extrasourceplots(self.variables[9], self.variables[10], self.eqcount, self.stcount,
                                                 self.variables[7], self.variables[5],
                                                 self.variables[4], self.variables[8])
            self.signalemit('Extra plots are saved', False)
        self.signalemit('All the selected modules are done.', True)
    def signalemit(self, text, style):
        self.signal.emit(text, style)

    def run(self):

        self.effectmethod()



class GeneralizedRun (QMainWindow):

    def __init__(self):

        QMainWindow.__init__(self)


    def initUi(self, mainwindow, variables):

        self.variables = variables


        # Setting main window changes
        mainwindow.setWindowTitle('Generalized Inversion')
        mainwindow.setGeometry(50, 50, 900, 700)

        # Setting layout for window after hitting run button
        mainwindow.centralWidget = QWidget(mainwindow)
        self.mainrunwidget = QWidget()

        self.mainrunlayout = QGridLayout()



# Widgets
        # Main text field creation
        self.runflow = QTextBrowser()
        self.runflow.isReadOnly()
        self.runflow.setText('The Program is now running:')
        self.runflow.moveCursor(QTextCursor.End)
        self.exitbtn = QPushButton('Exit')
        self.backbtn = QPushButton('Back')

# Adding widgets to layouts
        self.mainrunlayout.addWidget(self.runflow, 0, 0)
        self.mainrunlayout.addWidget(self.exitbtn, 1, 0)
        self.mainrunlayout.addWidget(self.backbtn, 2, 0)

        self.mainrunwidget.setLayout(self.mainrunlayout)
        mainwindow.centralWidget = self.mainrunwidget
        self.mainrunwidget.setWindowTitle('Generalized Inversion')
        self.mainrunwidget.setGeometry(50, 50, 500, 600)
        self.mainrunwidget.show()


        self.folder_making(variables[10])

        self.progress()

    def progress(self):
        self.process = EffectsCalculationThread(self.variables, self)
        self.process.start()
        # Connections
        self.process.signal.connect(self.updatetexteditor)
        self.backbtn.clicked.connect(self.backtohome)


    def updatetexteditor(self, text, style):

        header = QFont('Helvetica', 10)
        self.runflow.moveCursor(QTextCursor.End)
        self.runflow.setFont(header)


        if style:
            self.runflow.insertPlainText('\n\n' + text)
        else:
            self.runflow.insertPlainText('\n\n' + text + '\n' + '---------------------------------------' )



    def backtohome(self):
        self.mainrunwidget.close()

    def folder_making(self, outputfolder):

        if not os.path.exists(outputfolder):  # Making main reults folder
            os.makedirs(outputfolder)
        if not os.path.exists(outputfolder + 'Interpolating'):  # Making interpolating folder
            os.makedirs(outputfolder + 'Interpolating')
        if not os.path.exists(outputfolder + 'Geometrical-spreading'):  # Making Geometrical spreading folder
            os.makedirs(outputfolder + 'Geometrical-spreading')
        if not os.path.exists(outputfolder + 'matrices'):  # Making matrices folder
            os.makedirs(outputfolder + 'matrices')
        if not os.path.exists(outputfolder + 'svd'):  # Making SVD folder
            os.makedirs(outputfolder + 'svd')
        if not os.path.exists(outputfolder + 'svd/Covariance'):  # Making Covariance folder
            os.makedirs(outputfolder + 'svd/Covariance')
        if not os.path.exists(outputfolder + 'Path-effect'):  # Making Path results folder
            os.makedirs(outputfolder + 'Path-effect')
        if not os.path.exists(outputfolder + 'Siteamplifications'):  # Making site amplification results folder
            os.makedirs(outputfolder + 'Siteamplifications')
        if not os.path.exists(outputfolder + 'Source-effects'):  # maikng source effect results folder
            os.makedirs(outputfolder + 'Source-effects')
        if not os.path.exists(outputfolder + 'Source-effects/grid-search'):  # Grid search results folder
            os.makedirs(outputfolder + 'Source-effects/grid-search')
        if not os.path.exists(outputfolder + 'Source-effects/Extra-plots'):  # Grid search results folder
            os.makedirs(outputfolder + 'Source-effects/Extra-plots')