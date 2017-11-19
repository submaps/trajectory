import numpy as np
import logging


class StressCyl:
    def __init__(self, s_r, s_theta, s_z, t_r_theta, t_theta_z, t_r_z):
                       #s_r, s_theta, s_z, t_r_theta, t_theta_z, tr_z

        self.s_r = float(s_r)
        self.s_theta = float(s_theta)
        self.s_z = float(s_z)
        self.t_r_theta = float(t_r_theta)
        self.t_theta_z = float(t_theta_z)
        self.t_r_z = float(t_r_z)

    def getMtx33(self):
        m = [[self.s_r, self.t_r_theta, self.t_r_z],
             [self.t_r_theta, self.s_theta, self.t_theta_z],
             [self.t_r_z, self.t_theta_z, self.s_z]]
        # self.printStressCyl()

        return m
    #
    # def getGeneralStressMtx(self):
    #
    #     m = self.getMtx33()
    #     # s3, s2, s1 = np.linalg.eigvals(m)
    #     arr = np.linalg.eigvals(m)
    #     arr = sorted(arr, reverse=True)
    #     # s3 = (m[0][0]+m[1][1])/2+(((m[0][0]-m[1][1])**2)/2+m[0][1]**2)**0.5
    #     # s1 = (m[0][0] + m[1][1]) / 2 - (((m[0][0] - m[1][1]) ** 2) / 2 + m[0][1]**2) ** 0.5
    #     # s2=s1+1#!
    #     s3 = arr[0]
    #     s2 = arr[1]
    #     s1 = arr[2]
    #
    #     # if s3 > s2 and s2 > s1 or (s3 == 0 and s1 == 0):
    #     #     # logging.info("\t\tgeneral stress in atm:"+str(s3/100000)+" "+str(s2/100000)+" "+str(s1/100000))
    #     #     return arr
    #         # print "Ok!"
    #     # else:
    #         #logging.error("\t\t!s3>s2>s1: " + str(s3) + " " + str(s2) + " " + str(s1))
    #     return arr

    def getGeneralStress(self):

        m = self.getMtx33()
        arr = np.linalg.eigvals(m)
        arr = sorted(arr, reverse=True)

        # s1 = self.s_r
        #
        # s3 = (self.s_theta + self.s_z) / 2 + np.sqrt(
        #     (self.s_theta - self.s_z) * (self.s_theta - self.s_z) + 4 * self.t_theta_z * self.t_theta_z) / 2
        #
        # s2 = (self.s_theta + self.s_z) / 2 - np.sqrt(
        #     (self.s_theta - self.s_z) * (self.s_theta - self.s_z) + 4 * self.t_theta_z * self.t_theta_z) / 2

        # arr=[s1, s2, s3]
        #
        # arr = sorted(arr, reverse=True)

        s3 = arr[0]
        s2 = arr[1]
        s1 = arr[2]

        # if s3 >=s2 and s2 >= s1 or (s3 == 0 and s1 == 0):
        #     # logging.info("\t\tgeneral stress in atm:"+str(s3/100000)+" "+str(s2/100000)+" "+str(s1/100000))
        #     return arr
        #     # else:
            # logging.error("\t\t!s3>s2>s1: " + str(s3) + " " + str(s2) + " " + str(s1))
        return arr

    def getGeneralStressStr(self):
        gstress_str="\t"
        for stress in self.getGeneralStress():
            gstress_str+=str(stress)+"\t"
        return gstress_str

    def printStressCyl(self):
        # logging.info("s_r       =" + str(self.s_r))
        # logging.info("s_theta   =" + str(self.s_theta))
        # logging.info("s_z       =" + str(self.s_z))
        # logging.info("t_r_theta =" + str(self.t_r_theta))
        # logging.info("t_theta_z =" + str(self.t_theta_z))
        # logging.info("t_r_z     =" + str(self.t_r_z))
        print("_____________________________")
        print("s_r       =" + str(self.s_r))
        print("s_theta   =" + str(self.s_theta))
        print("s_z       =" + str(self.s_z))
        print("t_r_theta =" + str(self.t_r_theta))
        print("t_theta_z =" + str(self.t_theta_z))
        print("t_r_z     =" + str(self.t_r_z))
        print("_____________________________")

    def getMtx33Str(self):
        curStr = ""
        curStr += str(self.s_r) + "\t" + str(self.t_r_theta) + "\t" + str(self.t_r_z) + "\n"
        curStr += str(self.t_r_theta) + "\t" + str(self.s_theta) + "\t" + str(self.t_theta_z) + "\n"
        curStr += str(self.t_r_z) + "\t" + str(self.t_theta_z) + "\t" + str(self.s_z)
        return curStr

    def getTrace(self):
        return self.s_r+self.s_theta+self.s_z

    def get_gstress_trace(self):
        sigma3, sigma2, sigma1=self.getGeneralStress()
        return sigma3+sigma2+sigma1