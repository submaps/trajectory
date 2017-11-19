# -*- coding: utf-8 -*-
from trace_checker import*
import numpy as np
from config import*

def get_min_ro_in_section(j, i, a, sigma_H, sigma_h, sigma_v,
                          r_vector, theta_vector, gridShape,
                          Rw, z, nu, phi, C0):

    # l.z:3000.0,     h:500.0,     z:2636.619772,     point.z - point_prev.z = 0.0,
    # 62:        2636.619772, l:5;i_deg: 90.0;a_deg: 0.0
    # sigma_H: 369.211889539 atm sigma_h: 369.211889539 atm sigma_v: 628.658082188
    # pw:258.652399633 ro: 1000 Coulomb checking
    # 62.0 phi: 0.785398163397 C0: 25000000 \
    #     -1.0 x + y + -25000000
    # dist = 8457618.37032

    fcriteria = FailureCriteria(ro_max, j)

    for ro in np.arange(ro_min, ro_max, (ro_max-ro_min)/25):
        pw = ro * g * z
        cur_section_cyl = calculate_section_at_well(j, i, a, pw,
                                                  sigma_H, sigma_h, sigma_v,
                                                  r_vector, theta_vector, gridShape,
                                                  Rw, z, nu)
        # проверка критериями массива напряжений на стенке скважины
        fcriteria.checkSection(cur_section_cyl.stressAtWell, ro, j, C0, phi)
    # print(fcriteria.getRo())
    roList = fcriteria.getRoList()
    return roList

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