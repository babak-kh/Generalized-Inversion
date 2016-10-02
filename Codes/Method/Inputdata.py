import xlrd
import numpy as np
import matplotlib.pyplot as plt


class InputData(object):

    def earthquakedata(self):


        earthquakes = xlrd.open_workbook('/home/babak/PycharmProjects/Simulation/GeneralInversion/Inputs/'
                                    'earthquakes.xls')
        earthquakesheet = earthquakes.sheet_by_index(0)
        self.eqno = earthquakesheet.col_values(0)
        self.eqcordlat = earthquakesheet.col_values(7)
        self.eqcordlong = earthquakesheet.col_values(8)
        self.eqdepth = earthquakesheet.col_values(9)
        self.eqmag = earthquakesheet.col_values(12)
        self.earthquakedatalist = [self.eqno, self.eqcordlong, self.eqcordlat,
                                   self.eqdepth, self.eqmag]
        return self.earthquakedatalist


