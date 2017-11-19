# -*- coding: utf-8 -*-
import matplotlib.pyplot as pl
import numpy as np
def plot_failure_criteria(j, sigma3, sigma2, sigma1, phi, C0):
    if sigma3==0:
        return False

    fig = pl.figure()
    ax = fig.add_subplot(111,aspect='equal')

    X1 = range(int(-float(C0) / np.tan(phi)), int(sigma3), int(sigma3 / 20 + 1))
    Y1 = map(lambda x: np.tan(phi) * x + C0, X1)

    ax.plot(X1, Y1)
    # (x-x0)**2+(y-y0)**2=r**2
    # y=y0+(r**2-(x-x0)**2)**0.5
    if sigma3>sigma1:
        print "True"
    r = (sigma3 - sigma1) / 2.0

    print sigma3, " ", sigma2, " ", sigma1
    print r
    x0 = sigma1 + r
    y0 = 0
    X2 = range(int(sigma1), int(sigma3), int(sigma3 / 1000 + 1))
    Y2 = map(lambda x: y0 + (r ** 2 - (x - x0) ** 2) ** 0.5, X2)
    ax.plot(X2, Y2)

    # circle center
    ax.scatter(x0,y0,color="red")

    # left point
    ax.scatter(sigma1, 0)

    # right point
    ax.scatter(sigma3, 0,color="green")

    #highest point
    ax.scatter(sigma1 + r, r, color="yellow")

    # UCS
    ax.scatter(0, C0, color="gray")

    # Ax+By+C=0
    A = -np.tan(phi)
    B = 1
    C = -C0
    dist = np.abs(A * x0 + B * y0 + C)/(A*A+B*B)**0.5
    print "dist=", dist/100000
    print "r=", r/100000
    if (dist <= r):
        ans = "Failure occurs"
    else:
        ans = "Failure does not occurs"

    fig.text(
        0.4, 0.1,
        ans,
        horizontalalignment='left',
        fontsize=15,
        transform=ax.transAxes
    )
    pl.savefig("./output/failure/" + str(j) + ".png", dpi=200)
    # pl.show()
    pl.close()


def ro_plot(zMas, ro1_vector,
                  ro2_vector,
                  ro3_vector,
                  ro4_vector,
                  ro5_vector,
                  ro6_vector,
                  ro7_vector,
                  ro8_vector):
    fig=pl.figure()
    ax = fig.add_subplot(111, aspect='equal')
    # ax = fig.add_subplot(111, aspect='equal')
    pl.title(r'$\rho(z)$', size=14)
    pl.xlabel(r'$\rho$', size=14)
    pl.ylabel('z', size=14)

    ax.plot(ro1_vector, zMas, color="green", label="FailureCoulombPrev",  markevery=(0, 5))
    ax.plot(ro2_vector, zMas, color="red", label="FailureCoulomb")
    ax.plot(ro3_vector, zMas, color="blue", label="FailureTresk")
    ax.plot(ro4_vector, zMas, color="yellow", label="FailureMises")
    ax.plot(ro5_vector, zMas, color="black", label="CoulombMohr")
    ax.plot(ro6_vector, zMas, linestyle='--', label="DruckerPragerInnerCircle")
    ax.plot(ro7_vector, zMas, linestyle='-', label="DruckerPragerMiddleCircl")
    ax.plot(ro8_vector, zMas, label="DruckerPragerOuterCircle")

    # pl.legend(loc='upper left')
    pl.savefig("./output/ro_plot.png", dpi=200)
    pl.show()


def prepare_one_ax(ax, pl, theta_grid, r_grid, data, title, title_size, Rw, title_position_x, title_position_y):
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    pl.pcolormesh(theta_grid, r_grid, data)
    pl.colorbar(pad=0.1)
    # pl.clim(minval, maxval)
    ax.set_thetagrids(np.array([0, 90, 180, 270]), ['0', '90', '180', '270'])
    # ax.set_thetagrids(np.array([0, 90, 180, 270]), ['0', '90', '180', '270'],fontsize=8)
    # ax.set_rgrids(radii=[Rw * 1, Rw * 2], labels=['1', ' 2 '], angle=90, fontsize=15)
    # ax.set_rgrids(radii=[Rw * 1, Rw * 2], labels=['$R_w$', '$2R_w$'], angle=90, fontsize=8)
    ax.set_rgrids(radii=[Rw * 1, Rw * 2], angle=90, fontsize=8)
    # ax.set_rgrids(radii=[Rw * 1, Rw * 2], labels=['', ''], angle=90, fontsize=8)

    ax.grid(True, color='black', linestyle='-', linewidth=1, axis='y')
    pl.text(title_position_x, title_position_y,
            title,
            horizontalalignment='left',
            fontsize=title_size,
            transform=ax.transAxes)
    # ax.set_title(title)


def make_form_plot(theta_grid, r_grid, curSectionCyl, j):
    Rw = curSectionCyl.Rw
    i = curSectionCyl.i
    pw = curSectionCyl.pw
    z = curSectionCyl.z

    title1 = r'$\sigma_r$'
    title2 = r'$\sigma_\theta$'
    title3 = r'$\sigma_z$'
    title4 = r'$\tau_{r\theta}$'
    title5 = r'$\tau_{\theta z}$'
    title6 = r'$\tau_{rz}$'

    s_r_data = curSectionCyl.s_r
    s_theta_data = curSectionCyl.s_theta
    s_z_data = curSectionCyl.s_z
    t_r_theta_data = curSectionCyl.t_r_theta
    t_theta_z_data = curSectionCyl.t_theta_z
    t_r_z_data = curSectionCyl.t_r_z

    fig = pl.figure()
    title_size = 20

    title_position_x = -0.6
    title_position_y = 1

    ax1 = pl.subplot(331, polar=True)
    prepare_one_ax(ax1, pl, theta_grid, r_grid, s_r_data, title1, title_size, Rw, title_position_x, title_position_y)

    ax2 = pl.subplot(334, polar=True)
    prepare_one_ax(ax2, pl, theta_grid, r_grid, s_theta_data, title2, title_size, Rw, title_position_x, title_position_y)

    ax3 = pl.subplot(337, polar=True)
    prepare_one_ax(ax3, pl, theta_grid, r_grid, s_z_data, title3, title_size, Rw, title_position_x, title_position_y)

    ax4 = pl.subplot(333, polar=True)
    prepare_one_ax(ax4, pl, theta_grid, r_grid, t_r_theta_data, title4, title_size, Rw, title_position_x,
                   title_position_y)

    ax5 = pl.subplot(336, polar=True)
    prepare_one_ax(ax5, pl, theta_grid, r_grid, t_theta_z_data, title5, title_size, Rw, title_position_x,
                   title_position_y)

    ax6 = pl.subplot(339, polar=True)
    prepare_one_ax(ax6, pl, theta_grid, r_grid, t_r_z_data, title6, title_size, Rw, title_position_x, title_position_y)

    fig.text(
        0.45, 0.1,
        "i=" + str(round(i * 180 / np.pi, 2)) + "\npw=" + str(round(pw / 100000, 2)) + " atm\n" + str(z) + " m",
        horizontalalignment='left',
        fontsize=15,
        transform=ax1.transAxes
    )

    fig.set_label("Main graph")
    # pl.savefig("./tmp/111.png", dpi=200)
    # pl.legend()
    pl.show()
    pl.close()


def make_stress_plot(theta_grid, r_grid, cur_section_cyl, j):
    Rw = cur_section_cyl.Rw
    i = cur_section_cyl.i
    pw = cur_section_cyl.pw
    z = cur_section_cyl.z

    title1 = r'$\sigma_r$'
    title2 = r'$\sigma_\theta$'
    title3 = r'$\sigma_z$'
    title4 = r'$\tau_{r\theta}$'
    title5 = r'$\tau_{\theta z}$'
    title6 = r'$\tau_{rz}$'

    s_r_data = cur_section_cyl.s_r
    s_theta_data = cur_section_cyl.s_theta
    s_z_data = cur_section_cyl.s_z
    t_r_theta_data = cur_section_cyl.t_r_theta
    t_theta_z_data = cur_section_cyl.t_theta_z
    t_r_z_data = cur_section_cyl.t_r_z

    fig = pl.figure()
    title_size = 20

    ax = pl.subplot(331, polar=True)
    title_position_x = -0.6
    title_position_y = 1
    prepare_one_ax(ax, pl, theta_grid, r_grid, s_r_data, title1, title_size, Rw, title_position_x, title_position_y)

    ax = pl.subplot(334, polar=True)
    prepare_one_ax(ax, pl, theta_grid, r_grid, s_theta_data, title2, title_size, Rw, title_position_x, title_position_y)

    ax = pl.subplot(337, polar=True)
    prepare_one_ax(ax, pl, theta_grid, r_grid, s_z_data, title3, title_size, Rw, title_position_x, title_position_y)

    ax = pl.subplot(333, polar=True)
    prepare_one_ax(ax, pl, theta_grid, r_grid, t_r_theta_data, title4, title_size, Rw, title_position_x, title_position_y)

    ax = pl.subplot(336, polar=True)
    prepare_one_ax(ax, pl, theta_grid, r_grid, t_theta_z_data, title5, title_size, Rw, title_position_x, title_position_y)

    ax = pl.subplot(339, polar=True)
    prepare_one_ax(ax, pl, theta_grid, r_grid, t_r_z_data, title6, title_size, Rw, title_position_x, title_position_y)

    fig.text(
        0.45, 0.1,
        "i=" + str(round(i * 180 / np.pi, 2)) + "\npw=" + str(round(pw / 100000, 2)) + " atm\n" + str(z) + " m",
        horizontalalignment='left',
        fontsize=15,
        transform=ax.transAxes
    )

    fig.set_label("Main graph")
    pl.savefig('./output/' + str(j) + '.png', dpi=200)
    # pl.legend()
    # pl.show()
    pl.close()


def draw_points(raw_data):
    fig = pl.figure()
    ax = fig.add_subplot(111, aspect='equal')
    x_vector = list(raw_data.x)
    z_vector = list(-1 * raw_data.z)

    ax.plot(x_vector, z_vector, marker='^', linestyle='--')
    # ax.plot(x_vector, z_vector, marker='^', linestyle='--', markevery=(0, 5))
    text1 = pl.text(0.5, 0.44, u'Trajectory')
    ax.grid(True)  # линии вспомогательной сетки
    pl.savefig('./output/points.png', dpi=300)
    pl.show()
    pl.close()

# def plotFailureCriteria(j, sigma3, sigma2, sigma1, phi, C0):
#
#     if sigma3==0:
#         return False
#     fig = pl.figure()
#     ax = fig.add_subplot(111,aspect='equal')
#
#     X1 = range(int(-float(C0) / np.tan(phi)), int(sigma3), int(sigma3 / 20 + 1))
#     Y1 = map(lambda x: np.tan(phi) * x + C0, X1)
#
#     ax.plot(X1, Y1)
#     # (x-x0)**2+(y-y0)**2=r**2
#     # y=y0+(r**2-(x-x0)**2)**0.5
#     r = (sigma3 - sigma1) / 2.0
#
#     logging.info(" sigma3="+str(sigma3/100000)+
#                  " sigma2="+str(sigma2/100000)+
#                  " sigma1="+str(sigma1/100000)+
#                  " r="+str(r/100000))
#     x0 = sigma1 + r
#     y0 = 0
#     X2 = range(int(sigma1), int(sigma3), int(sigma3 / 1000 + 1))
#     Y2 = map(lambda x: y0 + (r ** 2 - (x - x0) ** 2) ** 0.5, X2)
#     ax.plot(X2, Y2)
#
#     #circle center
#     ax.scatter(x0,y0,color="red")
#
#     ax.scatter(sigma1,0)
#
#     #
#     ax.scatter(sigma3, 0, color="green")
#
#     #highest point
#     ax.scatter(sigma1 + r, r, color="yellow")
#
#     #UCS
#     ax.scatter(0, C0, color="gray")
#
#     pl.savefig("./output/failure" + str(j) + ".png", dpi=200)
#     # pl.show()
#     pl.close()


# l.z:1150.0 h:1000.0 z:960.0 l.h_k-(l.z - point.z)=810
# sigma_H: 1654.97035714 atm 	sigma_h: 1654.97035714	  atm sigma_v: 2817.9225	pw:12242880.0	ro: 1300.0

if __name__=='__main__':
    print 10 ** 5
    sigma3 = 2817 * 100000
    sigma2 = 2217 * 100000
    sigma1 = -1000 * 100000
    C0 = 20 * 1000000
    phi = 43 * np.pi / 180
    plot_failure_criteria(10, sigma3, sigma2, sigma1, phi, C0)
    # plotFailureCriteria(10, 0, 0, 0, phi, C0)

    ro1_vector = np.linspace(800, 900, 100)
    ro2_vector = np.linspace(900, 1100, 100)
    ro3_vector = np.linspace(1200, 1300, 100)
    zMas = np.linspace(0, -2636, 100)
    ro_plot(zMas, ro1_vector, ro2_vector, ro3_vector)




