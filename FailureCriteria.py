# -*- coding: utf-8 -*-
import numpy as np
import logging
from failureCriteriaGraph import plot_failure_criteria


class FailureCriteria:
    def __init__(self, ro_init, id):
        self.roFailureCoulombPrev =          ro_init
        self.roFailureCoulomb =              ro_init
        self.roFailureTresk =                ro_init
        self.roFailureMises =                ro_init
        self.roCoulombMohr =                 ro_init
        self.roDruckerPragerInnerCircle =    ro_init
        self.roDruckerPragerMiddleCircle =   ro_init
        self.roDruckerPragerOuterCircle =    ro_init
        self.ro_max = ro_init
        self.id = id

    def get_ro(self):
        roAsText = ""
        roAsText += "roFailureCoulombPrev           =" + str(self.roFailureCoulombPrev)+"\n"
        roAsText += "roFailureCoulomb               =" + str(self.roFailureCoulomb)+"\n"
        roAsText += "roFailureTresk                 =" + str(self.roFailureTresk)+"\n"
        roAsText += "roFailureMises                 =" + str(self.roFailureMises)+"\n"
        roAsText += "roCoulombMohr                  =" + str(self.roCoulombMohr)+"\n"
        roAsText += "roDruckerPragerInnerCircle     =" + str(self.roDruckerPragerInnerCircle)+"\n"
        roAsText += "roDruckerPragerMiddleCircle    =" + str(self.roDruckerPragerMiddleCircle)+"\n"
        roAsText += "roDruckerPragerOuterCircle     =" + str(self.roDruckerPragerOuterCircle)+"\n"
        return roAsText

    def getRoList(self):
        return [
            self.id,
            self.roFailureCoulombPrev,
            self.roFailureCoulomb,
            self.roFailureTresk,
            self.roFailureMises,
            self.roCoulombMohr,
            self.roDruckerPragerInnerCircle,
            self.roDruckerPragerMiddleCircle,
            self.roDruckerPragerOuterCircle
        ]

    def checkIsFailureCoulombPrev(self,j, sigma3, sigma2, sigma1, phi, C0):
        # sigma3>sigma2>sigma1
        # C0 прочность на одноосное сжатие
        if sigma3 == 0 and sigma2 == 0 and sigma1 == 0:
            # logging.info("all sigma equals zero")
            return False

        r = (sigma3 - sigma1) / 2
        x0 = sigma1 + r
        y0 = 0
        # Ax+By+C=0
        # logging.info("Coulomb checking " + str(j) + "phi: " + str(phi) + " C0: " + str(C0))
        A = -np.tan(phi)
        B = 1
        C = -C0
        # logging.info(str(A) + "x+" + str(B) + "y+" + str(C))
        # расстояние между центром окружоности и прямой
        dist = np.abs(A * x0 + B * y0 + C) / (A * A + B * B) ** 0.5
        # logging.info("dist=" + str(dist))
        if dist <= r:
            # logging.info("failure criteria plotting")
            # plotFailureCriteria(j, sigma3, sigma2, sigma1, phi, C0)
            return True
        else:
            return False

    def checkIsFailureCoulomb(self,j, sigma3, sigma2, sigma1, phi, C0):
        # sigma3>sigma2>sigma1
        # C0 прочность на одноосное сжатие
        # #phi угол внутреннего трения
        # phi = i*np.pi/180
        q = (1 + np.sin(phi)) / (1 - np.sin(phi))
        # ?
        ans = (sigma1 >= C0 + sigma3 * q)
        if ans:
            return True
        else:
            return False

    def checkIsFailureTresk(self,j, sigma3, sigma2, sigma1, phi, C0):
        # sigma3>sigma2>sigma1
        # C_0 прочность на одноосное сжатие
        # phi угол внутреннего трения

        phi = 0
        # ?
        q = (1 + np.sin(phi)) / (1 - np.sin(phi))
        ans = (sigma1 >= C0 + sigma3 * q)
        if ans:
            return True
        else:
            return False

    def checkIsFailureMises(self,j, sigma3, sigma2, sigma1, phi, C0):
        # sigma3>sigma2>sigma1
        # sigma_0 прочность на одноосное сжатие
        # #phi угол внутреннего трения
        J2 = 1 / 6 * ((sigma1 - sigma2) ** 2 + (sigma1 - sigma3) ** 2 + (sigma2 - sigma3) ** 2)
        ans = (J2 ** 0.5 - C0 / 3 > 0)
        # ans=(J2**0.5-C0/3>0)
        # условие разрушения
        if ans:
            return True
        else:
            return False

    def checkCoulombMohr(self,j, sigma3, sigma2, sigma1, phi, C0):
        # return True if the failure occurs
        A = (1 + np.sin(phi)) / (1 - np.sin(phi))
        B = 2 * C0 * np.cos(phi) / (1 - np.sin(phi))
        f = A * sigma3 + B - sigma1
        # plotFailureCriteria(j, sigma3, sigma2, sigma1, phi, C0)
        if (f < 0):
            return True
        return False

    def checkDruckerPragerInnerCircle(self,j, sigma3, sigma2, sigma1, phi, C0):
        # return True if the failure occurs
        A = 6 ** 0.5 * np.sin(phi) / (9 + 3 * np.sin(phi) ** 2) ** 0.5
        B = 6 ** 0.5 * C0 * np.sin(phi) / (9 + 3 * np.sin(phi) ** 2) ** 0.5
        J1 = sigma1 + sigma2 + sigma3
        J2 = sigma2 * sigma3 + sigma1 * sigma3 + sigma1 * sigma2
        f = A * J1 + B - J2 ** 2
        if (f < 0):
            return True
        return False

    def checkDruckerPragerMiddleCircle(self,j, sigma3, sigma2, sigma1, phi, C0):
        # return True if the failure occurs
        A = 2 * 2 ** 2 * np.sin(phi) / (3 + np.sin(phi))
        B = 2 * 2 ** 2 * np.cos(phi) / (3 + np.sin(phi))
        J1 = sigma1 + sigma2 + sigma3
        J2 = sigma2 * sigma3 + sigma1 * sigma3 + sigma1 * sigma2
        f = A * J1 + B - J2 ** 2
        if (f < 0):
            return True
        return False

    def checkDruckerPragerOuterCircle(self,j, sigma3, sigma2, sigma1, phi, C0):
        # return True if the failure occurs
        A = 2 * 2 ** 2 * np.sin(phi) / (3 - np.sin(phi))
        B = 2 * 2 ** 2 * np.cos(phi) / (3 - np.sin(phi))
        J1 = sigma1 + sigma2 + sigma3
        J2 = sigma2 * sigma3 + sigma1 * sigma3 + sigma1 * sigma2
        f = A * J1 + B - J2 ** 2
        if (f < 0):
            return True
        return False

    def checkDruckerPragerPlainStrainState(self,j, sigma3, sigma2, sigma1, phi, C0):
        # return True if the failure occurs
        A = 2 * np.sin(phi) / (9 - 3 * np.sin(phi)) ** 0.5
        B = 6 * C0 * np.sin(phi) / (9 - 3 * np.sin(phi)) ** 0.5
        J1 = sigma1 + sigma2 + sigma3
        J2 = sigma2 * sigma3 + sigma1 * sigma3 + sigma1 * sigma2
        f = A * J1 + B - J2 ** 2
        if f < 0:
            return True
        return False

    def checkDruckerPragerTriaxialCompression(self,j, sigma3, sigma2, sigma1, phi, C0):
        # return True if the failure occurs
        A = np.tan(phi) / (9 + 12 * np.tan(phi) ** 2) ** 0.5
        B = 3 * C0 / (9 + 12 * np.tan(phi) ** 2) ** 0.5
        J1 = sigma1 + sigma2 + sigma3
        J2 = sigma2 * sigma3 + sigma1 * sigma3 + sigma1 * sigma2
        f = A * J1 + B - J2 ** 2
        if (f < 0):
            return True
        return False

    def checkSection(self, stressAtWell, ro, j, C0, phi):
        # главные напряжения в этой точке
        # sigma3>sigma2>sigma1
        # phi = 45.0 * np.pi / 180
        # C0 = 25 * 1000000

        # по умолчанию все критерии не выполняются
        isCheckIsFailureCoulombPrev = False
        isCheckIsFailureCoulomb = False
        isCheckIsFailureTresk = False
        isCheckIsFailureMises = False
        isCheckCoulombMohr = False
        isCheckDruckerPragerInnerCircle = False
        isCheckDruckerPragerMiddleCircle = False
        isCheckDruckerPragerOuterCircle = False
        isCheckDruckerPragerPlainStrainState = False

        # print("checkSection started\n")

        # print stressAtWell.__sizeof__()

        # print("===========================================")

        for curNumber, currentStressCyl in enumerate(stressAtWell):
            # print(str(curNumber) + ":")
            # print(currentStressCyl.getMtx33Str())
            sigma3, sigma2, sigma1 = currentStressCyl.getGeneralStress()

            # print(sigma3>sigma2, sigma2>sigma1, sigma3>sigma1)
            # print("\ts3=" + str(sigma3) + "\ts2=" + str(sigma2) + "\ts1=" + str(sigma1))

            if isCheckIsFailureCoulombPrev != True:
                isCheckIsFailureCoulombPrev = self.checkIsFailureCoulombPrev(curNumber, sigma3, sigma2, sigma1, phi, C0)

            if isCheckIsFailureCoulomb != True:
                isCheckIsFailureCoulomb = self.checkIsFailureCoulomb(curNumber, sigma3, sigma2, sigma1, phi, C0)

            if isCheckIsFailureTresk != True:
                isCheckIsFailureTresk = self.checkIsFailureTresk(curNumber, sigma3, sigma2, sigma1, phi, C0)

            if isCheckIsFailureMises != True:
                isCheckIsFailureMises = self.checkIsFailureMises(curNumber, sigma3, sigma2, sigma1, phi, C0)

            if isCheckCoulombMohr != True:
                isCheckCoulombMohr = self.checkCoulombMohr(curNumber, sigma3, sigma2, sigma1, phi, C0)

            if isCheckDruckerPragerInnerCircle != True:
                isCheckDruckerPragerInnerCircle = self.checkDruckerPragerInnerCircle(curNumber, sigma3, sigma2, sigma1,
                                                                                     phi,
                                                                                     C0)
            if isCheckDruckerPragerMiddleCircle != True:
                isCheckDruckerPragerMiddleCircle = self.checkDruckerPragerMiddleCircle(curNumber, sigma3, sigma2,
                                                                                       sigma1,
                                                                                       phi, C0)
            if isCheckDruckerPragerOuterCircle != True:
                isCheckDruckerPragerOuterCircle = self.checkDruckerPragerOuterCircle(curNumber, sigma3, sigma2, sigma1,
                                                                                     phi, C0)
            if isCheckDruckerPragerPlainStrainState != True:
                isCheckDruckerPragerPlainStrainState = self.checkDruckerPragerPlainStrainState(curNumber, sigma3,
                                                                                               sigma2,
                                                                                               sigma1,
                                                                                               phi, C0)
        if isCheckIsFailureCoulombPrev == False:
            if ro < self.roFailureCoulombPrev:
                self.roFailureCoulombPrev = ro

        if isCheckIsFailureCoulomb == False:
            if ro < self.roFailureCoulomb:
                self.roFailureCoulomb = ro

        if isCheckIsFailureTresk == False:
            if ro < self.roFailureTresk:
                self.roFailureTresk = ro

        if isCheckIsFailureMises == False:
            if ro < self.roFailureMises:
                self.roFailureMises = ro

        if isCheckCoulombMohr == False:
            if ro < self.roCoulombMohr:
                self.roCoulombMohr = ro

        if isCheckDruckerPragerInnerCircle == False:
            if ro < self.roDruckerPragerInnerCircle:
                self.roDruckerPragerInnerCircle = ro

        if isCheckDruckerPragerMiddleCircle == False:
            if ro < self.roDruckerPragerMiddleCircle:
                self.roDruckerPragerMiddleCircle = ro

        if isCheckDruckerPragerOuterCircle == False:
            if ro < self.roDruckerPragerOuterCircle:
                self.roDruckerPragerOuterCircle = ro

        # print('checkIsFailureCoulombPrev         ' + str(isCheckIsFailureCoulombPrev))
        # print('checkIsFailureCoulomb             ' + str(isCheckIsFailureCoulomb))
        # print('checkIsFailureTresk               ' + str(isCheckIsFailureTresk))
        # print('checkIsFailureMises               ' + str(isCheckIsFailureMises))
        # print('checkCoulombMohr                  ' + str(isCheckCoulombMohr))
        # print('checkDruckerPragerInnerCircle     ' + str(isCheckDruckerPragerInnerCircle))
        # print('checkDruckerPragerMiddleCircle     ' + str(isCheckDruckerPragerMiddleCircle))
        # print('checkDruckerPragerOuterCircle     ' + str(isCheckDruckerPragerOuterCircle))
        # print("---------------------------------------------------------------------------")

        # logging.info('checkIsFailureCoulombPrev         ' + str(isCheckIsFailureCoulombPrev))
        # logging.info('checkIsFailureCoulomb             ' + str(isCheckIsFailureCoulomb))
        # logging.info('checkIsFailureTresk               ' + str(isCheckIsFailureTresk))
        # logging.info('checkIsFailureMises               ' + str(isCheckIsFailureMises))
        # logging.info('checkCoulombMohr                  ' + str(isCheckCoulombMohr))
        # logging.info('checkDruckerPragerInnerCircle     ' + str(isCheckDruckerPragerInnerCircle))
        # logging.info('checkDruckerPragerMiddleCircle     ' + str(isCheckDruckerPragerMiddleCircle))
        # logging.info('checkDruckerPragerOuterCircle     ' + str(isCheckDruckerPragerOuterCircle))
