# -*- coding: utf-8 -*-
from trace_checker import*
import numpy as np
from config import *
from failure_criteria_graph import *
from trace_checker import calculate_section_at_well


def get_min_ro_in_section(j, i, a, sigma_H, sigma_h, sigma_v,
                          r_vector, theta_vector, gridShape,
                          Rw, z, nu, phi, C0):
    # объект, хранящий в себе все критерии
    # fcriteria = FailureCriteria(RO_MAX, j)
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
        fcriteria_max_ro_list = [RO_MIN]*len(fcriteria_ans)
        # 1300 800 1300 800 800 1300 1300 1300 1300 1300 1300
        for id, criterion in enumerate(fcriteria.criteria_list):
            # if not fcriteria_ans[criterion_id]:
            #     fcriteria_ro_list[criterion_id] = ro
            #если обрушения возникает то минимум такой
            if not fcriteria_ans[id]:
                fcriteria_ro_list[id] = ro
                # if ro > fcriteria_ro_list[criterion_id]:
                #     fcriteria_ro_list[criterion_id] = ro

            # if not fcriteria_ans[criterion_id]:
            #     if fcriteria_ro_list[criterion_id] == RO_MIN:
            #         fcriteria_ro_list[criterion_id] = ro
            #     else:
            #         if ro < fcriteria_ro_list[criterion_id]:
            #             fcriteria_ro_list[criterion_id] = ro

        # fcriteria.check_section(cur_section_cyl.stress_at_well, ro, j, phi, C0)
    # print(fcriteria.getRo())
    # roList = fcriteria.getRoList()
    return fcriteria_ro_list

# def get_min_ro_in_section_test(j, i, a, sigma_H, sigma_h, sigma_v,
#                               r_vector, theta_vector, gridShape,
#                               Rw, z, nu, phi, C0):
#         # объект, хранящий в себе все критерии
#         fcriteria = FailureCriteria(RO_MAX, j)
#         from trace_checker import calculate_section_at_well
#         for ro in np.arange(RO_MIN, RO_MAX, (RO_MAX - RO_MIN) / 25):
#             pw = ro * g * z
#             cur_section_cyl = calculate_section_at_well(j, i, a, pw,
#                                                         sigma_H, sigma_h, sigma_v,
#                                                         r_vector, theta_vector, gridShape,
#                                                         Rw, z, nu)
#
#             # проверка критериями массива напряжений на стенке скважины
#             fcriteria.check_section(cur_section_cyl.stress_at_well, ro, j, phi, C0)
#         # print(fcriteria.getRo())
#         roList = fcriteria.getRoList()
#         return roList

# def get_min_ro_in_section(j, i, a, sigma_H, sigma_h, sigma_v,
#                               r_vector, theta_vector, gridShape,
#                               Rw, z, nu, phi, C0):
#
#
#     # l.z:3000.0,     h:500.0,     z:2636.619772,     point.z - point_prev.z = 0.0,
#     # 62:        2636.619772, l:5;i_deg: 90.0;a_deg: 0.0
#     # sigma_H: 369.211889539 atm sigma_h: 369.211889539 atm sigma_v: 628.658082188
#     # pw:258.652399633 ro: 1000 Coulomb checking
#     # 62.0 phi: 0.785398163397 C0: 25000000 \
#     #     -1.0 x + y + -25000000
#     # dist = 8457618.37032
#
#     fcriteria = FailureCriteria(RO_MAX, j)
#
#     for ro in np.arange(RO_MIN, RO_MAX, (RO_MAX - RO_MIN) / 25):
#         pw = ro * g * z
#         cur_section_cyl = calculate_section_at_well(j, i, a, pw,
#                                                     sigma_H, sigma_h, sigma_v,
#                                                     r_vector, theta_vector, gridShape,
#                                                     Rw, z, nu)
#         # проверка критериями массива напряжений на стенке скважины
#         fcriteria.check_section(cur_section_cyl.stress_at_well, ro, j, phi, C0)
#     # print(fcriteria.getRo())
#     roList = fcriteria.getRoList()
#     return roList

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