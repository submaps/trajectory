# -*- coding: utf-8 -*-
import random

import numpy as np
import logging
from failure_criteria_graph import plot_failure_criteria
import stress_cyl
from config import *

class Criterion:
    def __init__(self, ro_init, name):
        self.ro = ro_init
        self.name = name

    def checkThisCriterion(self, j, sigma3, sigma2, sigma1, phi, C0):
        if self.name == "CoulombMohr":       return self.checkFailureCoulombPrev(j, sigma3, sigma2, sigma1, phi, C0)
        if self.name == "FailureCoulomb":           return self.checkFailureCoulomb(j, sigma3, sigma2, sigma1, phi, C0)
        if self.name == "Tresk":             return self.checkFailureTresk(j, sigma3, sigma2, sigma1, phi, C0)
        if self.name == "FailureMises":             return self.checkFailureMises(j, sigma3, sigma2, sigma1, phi, C0)
        if self.name == "CoulombMohrP":              return self.checkCoulombMohr(j, sigma3, sigma2, sigma1, phi, C0)
        if self.name == "DruckerPragerInnerCircle": return self.checkDruckerPragerInnerCircle(j, sigma3, sigma2, sigma1, phi, C0)
        if self.name == "DruckerPragerMiddleCircle":return self.checkDruckerPragerMiddleCircle(j, sigma3, sigma2, sigma1, phi, C0)
        if self.name == "DruckerPragerOuterCircle": return self.checkDruckerPragerOuterCircle(j, sigma3, sigma2, sigma1, phi, C0)

        if self.name == "DruckerPragerPlainStrain":return self.checkDruckerPragerPlainStrainState(j, sigma3, sigma2, sigma1, phi, C0)
        if self.name == "DruckerPragerTriaxial":  return self.checkDruckerPragerTriaxialCompression(j, sigma3, sigma2, sigma1, phi, C0)
        if self.name == "TestCriteria":  return self.testCriteria(j, sigma3, sigma2, sigma1, phi, C0)
        if self.name == "MogiCoulomb":  return self.checkMogiCoulomb(j, sigma3, sigma2, sigma1, phi, C0)
        if self.name == "ModifiedLade":  return self.checkModifiedLade(j, sigma3, sigma2, sigma1, phi, C0)

        return False

    def check_cur_ro(self, j, ro, sigma3, sigma2, sigma1, phi, C0):
        if not self.checkThisCriterion(j, ro, sigma3, sigma2, sigma1, phi, C0):
            if ro < self.ro:
                self.ro = ro
        # else:
            # logging.info("\t\t\tcriterion reached ro="+str(self.ro)+"\t"+self.name)

    def checkFailureCoulombPrev(self, j, sigma3, sigma2, sigma1, phi, C0):
        # sigma3>sigma2>sigma1
        # C0 прочность на одноосное сжатие
        if sigma3 == 0 and sigma2 == 0 and sigma1 == 0:
            # logging.info("all sigma equals zero")
            return True

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

        # plot_failure_criteria(j, sigma3, sigma2, sigma1, phi, C0)

        if r>=dist:
            # logging.info("failure criteria plotting")
            return False
        else:
            return True

    def checkFailureCoulomb(self, j, s1, s2, s3, phi, C0):
        # !!! s1>s2>s3
        # C0 прочность на одноосное сжатие
        # #phi угол внутреннего трения
        # phi = i*np.pi/180
        q = (1 + np.sin(phi)) / (1 - np.sin(phi))
        # ?
        ans = (s3 >= (C0 + s1 * q))
        if ans:
            return False
        else:
            return True

    def checkFailureTresk(self, j, s1, s2, s3, phi, C0):
        # sigma3>sigma2>sigma1
        # C_0 прочность на одноосное сжатие
        # phi угол внутреннего трения

        ans = s1-s3 - C0
        if ans>0:
            return False
        else:
            return True

    def testCriteria(self, j, sigma3, sigma2, sigma1, phi, C0):
        # sigma3>sigma2>sigma1
        # C_0 прочность на одноосное сжатие
        # phi угол внутреннего трения

        ans = random.random()
        if ans:
            return False
        else:
            return True

    def checkFailureMises(self, j, s1, s2, s3, phi, C0):

        # sigma_0 прочность на одноосное сжатие
        # #phi угол внутреннего трения
        # !!! s1>s2>s3
        J2 = 1 / 6 * ((s1 - s2) ** 2 + (s1 - s3) ** 2 + (s2 - s3) ** 2)
        ans = (J2 ** 0.5 - C0 / 3 > 0)
        # ans=(J2**0.5-C0/3>0)
        # условие разрушения
        if ans:
            return False
        else:
            return True

    def checkCoulombMohr(self, j, sigma3, sigma2, sigma1, phi, C0):
        # return False if the failure occurs
        A = (1 + np.sin(phi)) / (1 - np.sin(phi))
        B = 2 * C0 * np.cos(phi) / (1 - np.sin(phi))
        f = A * sigma3 + B - sigma1
        # plotFailureCriteria(j, sigma3, sigma2, sigma1, phi, C0)
        if f < 0:
            return False
        return True

    def checkDruckerPragerInnerCircle(self, j, sigma3, sigma2, sigma1, phi, C0):
        # return False if the failure occurs
        A = 6 ** 0.5 * np.sin(phi) / (9 + 3 * np.sin(phi) ** 2) ** 0.5
        B = 6 ** 0.5 * C0 * np.sin(phi) / (9 + 3 * np.sin(phi) ** 2) ** 0.5
        J1 = sigma1 + sigma2 + sigma3
        # J2 = ((sigma1-sigma2)**2+(sigma2-sigma3)**2+(sigma3-sigma1)**2)/6
        J2 = ((sigma1-sigma2)**2+(sigma2-sigma3)**2+(sigma1-sigma3)**2)/6
        f = A * J1 + B - J2 ** 0.5
        if f < 0:
            return False
        return True

    def checkDruckerPragerMiddleCircle(self, j, sigma3, sigma2, sigma1, phi, C0):
        # return False if the failure occurs
        A = 2 * 2 ** 0.5 * np.sin(phi) / (3 + np.sin(phi))
        B = 2 * 2 ** 0.5 * C0 * np.cos(phi) / (3 + np.sin(phi))
        J1 = sigma1 + sigma2 + sigma3
        J2 = ((sigma1 - sigma2) ** 2 + (sigma2 - sigma3) ** 2 + (sigma3 - sigma1) ** 2) / 6
        f = A * J1 + B - J2 ** 0.5
        if f < 0:
            return False
        return True

    def checkDruckerPragerOuterCircle(self, j, sigma3, sigma2, sigma1, phi, C0):
        # return False if the failure occurs
        A = 2 * 2 ** 0.5 * np.sin(phi) / (3 - np.sin(phi))
        B = 2 * 2 ** 0.5 * C0* np.cos(phi) / (3 - np.sin(phi))
        J1 = sigma1 + sigma2 + sigma3
        J2 = ((sigma1 - sigma2) ** 2 + (sigma2 - sigma3) ** 2 + (sigma3 - sigma1) ** 2) / 6
        f = A * J1 + B - J2 ** 0.5
        if (f < 0):
            return False
        return True

    def checkDruckerPragerPlainStrainState(self, j, sigma3, sigma2, sigma1, phi, C0):
        # return False if the failure occurs
        A = 2 * np.sin(phi) / (9 - 3 * np.sin(phi)) ** 0.5
        B = 6 * C0 * np.sin(phi) / (9 - 3 * np.sin(phi)) ** 0.5
        J1 = sigma1 + sigma2 + sigma3
        J2 = ((sigma1 - sigma2) ** 2 + (sigma2 - sigma3) ** 2 + (sigma3 - sigma1) ** 2) / 6
        f = A * J1 + B - J2 ** 0.5
        if f < 0:
            return False
        return True

    def checkDruckerPragerTriaxialCompression(self, j, sigma3, sigma2, sigma1, phi, C0):
        # return False if the failure occurs
        A = np.tan(phi) / (9 + 12 * np.tan(phi) ** 2) ** 0.5
        B = 3 * C0 / (9 + 12 * np.tan(phi) ** 2) ** 0.5
        J1 = sigma1 + sigma2 + sigma3
        J2 = ((sigma1 - sigma2) ** 2 + (sigma2 - sigma3) ** 2 + (sigma3 - sigma1) ** 2) / 6
        f = A * J1 + B - J2 ** 0.5
        if f < 0:
            return False
        return True

    def checkMogiCoulomb(self, j, s1, s2, s3, phi, C0):
        # !!! s1>s2>s3
        q = (1+np.sin(phi))/(1-np.sin(phi))
        a = 2*2**0.5/3 * C0/(q+1)
        b = 2 * 2 ** 0.5 / 3 * (q-1) / (q + 1)
        sigma_m2=(s1+s2)/2

        f = a*3/(2*2**0.5)+b*3/(2*2**0.5)*sigma_m2-(s1-3)/2

        # t_max = (a * 3 / (2 * 2 ** 0.5)) + (b * 3 / (2 * 2 ** 0.5)) * sigma_m2

        # t_oct=1/3*((sigma1-sigma2)**2+(sigma1-sigma3)**2+(sigma2-sigma3)**2)**0.5
        # t_oct=a+b*sigma_m2

        if f < 0:
            return False
        return True
    def checkModifiedLade(self, j, sigma3, sigma2, sigma1, phi, C0):

        S = C0 / np.tan(phi)
        ng=4*np.tan(phi)**2*(9-7*np.sin(phi))/(1-np.sin(phi))
        I1__= (sigma1+S)+(sigma2+S)+(sigma3+S)
        I3__= (sigma1+S)*(sigma2+S)*(sigma3+S)
        f=-I1__**3/I3__+27+ng

        if f < 0:
            return False
        return True



class FailureCriteria:
    def __init__(self, ro_init, id):
        self.id = id
        self.criteria_list = []

        self.criteria_names = [
            "CoulombMohr",#0
            # "FailureCoulomb",      #1
            "Tresk",      #2
            # "FailureMises",      #3
            # "CoulombMohrP",         #4
            # "DruckerPragerInnerCircle", #5
            # "DruckerPragerMiddleCircle", #6
            # "DruckerPragerOuterCircle", #7
            "DruckerPragerPlainStrain", #8
            "DruckerPragerTriaxial", #9
            "MogiCoulomb", #10
            "ModifiedLade" #11
            # "TestCriteria" #12
        ]

        for name in self.criteria_names:
            self.criteria_list.append(Criterion(ro_init, name))

        # массив ответов по всему сечению
        self.ans_list=[True]*len(self.criteria_list)

    # возвращает лист с решениями по каждому критерию в сечении
    # False = Failure
    def check_section(self, stress_at_well, phi, C0):
        # перед проверкой сечения все ответы True
        self.ans_list = [True]*len(self.criteria_list)
        for stress_id, stress_cyl in enumerate(stress_at_well):
            # logging.info(str(stress_id) + ": " + stress_cyl.getGeneralStressStr())
            # logging.info("\ttrace comparison: " + str(stress_cyl.getTrace() - stress_cyl.get_gstress_trace()))
            # logging.info("-------------------------------")
            # проверяем по каждому критерию текующую точку и прибавляем к листу ответов по сечению
            self.check_gstress_at_section(stress_cyl.getGeneralStress(),
                                                    phi,
                                                    C0)
            # sigma3, sigma2, sigma1 = stress_cyl.getGeneralStress()
            # plot_failure_criteria(stress_id, sigma3, sigma2, sigma1, phi, C0)
            # logging.info("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        # logging.info("criterion_ans_at_section:")
        # logging.info(self.ans_list)
        return self.ans_list

    # в заданной точке
    def check_gstress_at_point(self, gstress, phi, C0):
        sigma3, sigma2, sigma1 = gstress

        # в начале ответ по всем критериям True
        # как только по одному критерию срабатывает ответ False/Fail
        # решение по нему False
        fcriteria_ans = [True for c in self.criteria_list]
        for id, criterion in enumerate(self.criteria_list):
            fcriteria_ans[id] = fcriteria_ans[id] and criterion.checkThisCriterion(id, sigma3, sigma2, sigma1, phi, C0)
        # в результате получаем лист ответов по каждому критерию
        # по данном тензору главных напряжений в одной точке
        return fcriteria_ans

    def check_gstress_at_section(self, gstress, phi, C0):
        sigma3, sigma2, sigma1 = gstress

        # в начале ответ по всем критериям True
        # как только по одному критерию срабатывает ответ False/Fail
        # решение по нему False

        for id, criterion in enumerate(self.criteria_list):
            self.ans_list[id] = self.ans_list[id] and criterion.checkThisCriterion(id, sigma3, sigma2, sigma1, phi, C0)
        # в результате получаем лист ответов по каждому критерию
        # по данном тензору главных напряжений в одной точке
        return self.ans_list

