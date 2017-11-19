import logging
import numpy as np


class StressDec:
    def __init__(self, s_x, s_y, s_z, t_y_z, t_x_z, t_x_y):
                      #s_x, s_y, s_z, t_y_z, t_x_z, t_x_y
        self.s_x = float(s_x)
        self.s_y = float(s_y)
        self.s_z = float(s_z)
        self.t_y_z = float(t_y_z)
        self.t_x_z = float(t_x_z)
        self.t_x_y = float(t_x_y)

    def printStressDec(self):
        logging.info("============Stress at well===========")
        logging.info("\t\t\ts_x   =" + str(self.s_x))
        logging.info("\t\t\ts_y   =" + str(self.s_y))
        logging.info("\t\t\ts_z   =" + str(self.s_z))
        logging.info("\t\t\tt_y_z =" + str(self.t_y_z))
        logging.info("\t\t\tt_x_z =" + str(self.t_x_z))
        logging.info("\t\t\tt_x_y =" + str(self.t_x_y))
        logging.info("=======================")

    def getMtx33(self):
        m = [[self.s_x, self.t_x_y, self.t_x_z],
             [self.t_x_y, self.s_y, self.t_y_z],
             [self.t_x_z, self.t_y_z, self.s_z]]
        return m

    # def getGeneralStress(self):
    #     m = self.getMtx33()
    #     # s3, s2, s1 = np.linalg.eigvals(m)
    #     arr = np.linalg.eigvals(m)
    #     arr = sorted(arr, reverse=True)
    #     s3 = arr[0]
    #     s2 = arr[1]
    #     s1 = arr[2]
    #     if s3 > s2 and s2 > s1 or (s3 == 0 and s2 == 0 and s1 == 0):
    #         # logging.info(self.getMtx33())
    #         return arr
    #         # print "Ok!"
    #     # else:
    #     #     logging.error("!s3>s2>s1: " + str(s3) + " " + str(s2) + " " + str(s1))
    #     return arr

    def getTrace(self):
        return self.s_x+self.s_y+self.s_z
