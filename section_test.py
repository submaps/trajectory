# -*- coding: utf-8 -*-
from trace_checker import*
import numpy as np
from config import *
from failure_criteria_graph import *
from trace_checker import calculate_section_at_well

# возвращает минимальную плотность по каждому критерию
def get_min_ro_in_section(j, i, a, sigma_H, sigma_h, sigma_v,
                          r_vector, theta_vector, gridShape,
                          Rw, z, nu, phi, C0):
    # объект, хранящий в себе все критерии
    # fcriteria = FailureCriteria(RO_MAX, j)
    # j идентификатор
    fcriteria = FailureCriteria(RO_MIN, j)
    fcriteria_ro_list = [RO_MIN]*len(fcriteria.criteria_list)
    for ro in np.arange(RO_MIN, RO_MAX, (RO_MAX-RO_MIN)/25):
        pw = ro * g * z
        stress_at_well_list = calculate_section_at_well(j, i, a, pw,
                                                  sigma_H, sigma_h, sigma_v,
                                                  r_vector, theta_vector, gridShape,
                                                  Rw, z, nu)
        # проверка критериями массива напряжений на стенке скважины
        #[True True True False False True False False]
        #False~Failure
        fcriteria_ans = fcriteria.check_section(stress_at_well_list, phi, C0)
        logging.info( "\t\t"+str(ro)+":\t"+str(fcriteria_ans))
        for id, criterion in enumerate(fcriteria.criteria_list):
            #если обрушение возникает то минимум равен текущему
            if not fcriteria_ans[id]:
                fcriteria_ro_list[id] = ro
    return fcriteria_ro_list

if __name__=='__main__':
    sigma_H = 369.211889539 * 100000
    sigma_h = 369.211889539 * 100000
    sigma_v = -628.658082188 * 100000

    nu = 0.34

    j = 12345
    Rw = 0.01
    i = 45 * np.pi / 180
    a = 0
    z = 2636

    theta_vector = np.linspace(0, 2 * np.pi, 100)
    r_vector = np.linspace(Rw, Rw * 2, 25)
    r_grid, theta_grid = np.meshgrid(r_vector, theta_vector)
    gridShape = r_grid.shape

    get_min_ro_in_section(j, i, a, sigma_H, sigma_h, sigma_v,
                          r_vector, theta_vector, gridShape,
                          Rw, z, nu)

    # make_form_plot(
    #     theta_grid,
    #     r_grid,
    #     curSectionCyl,
    #     j
    # )

    # for rc, r in enumerate(r_vector):
    #     print(str(rc)+": "+ str(r))
    #     for thetac, theta in enumerate(theta_vector):
    #         pass

    roMas=[]
    for a in range(0,10):
        roMas.append([z, 800, 100000, 800, 100000])