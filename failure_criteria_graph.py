# -*- coding: utf-8 -*-
import matplotlib.pyplot as pl
import numpy as np
import logging
from config import *
import random


def plot_failure_criteria(id, sigma3, sigma2, sigma1, phi, C0):
    # sigma3>sigma2>sigma1
    if sigma3 == 0:
        return False
    fig = pl.figure()
    ax = fig.add_subplot(111,aspect='equal')
    # ax = fig.add_subplot(111)

    # диапазон X1 от точки пересечения прямой с осью х, точка (T0,0)
    X1 = range(int(-float(C0) / np.tan(phi)), int(sigma3), int(sigma3 / 20 + 1))
    Y1 = map(lambda x: np.tan(phi) * x + C0, X1)

    ax.plot(X1, Y1, color='red')
    # (x-x0)**2+(y-y0)**2=r**2
    # y=y0+(r**2-(x-x0)**2)**0.5
    if sigma3 > sigma1:
        logging.info("sigma3>sigma1: True")
    else:
        logging.info("sigma3>sigma1: False")
    r = (sigma3 - sigma1) / 2.0

    # print sigma3, " ", sigma2, " ", sigma1
    # print "\tr: ", r
    x0 = sigma1 + r
    y0 = 0
    X2 = range(int(sigma1), int(sigma3), int(sigma3 / 1000 + 1))
    # Y2 = map(lambda x: y0 + (r ** 2 - (x - x0) ** 2) ** 0.5, X2)

    Y2 = [y0+(r ** 2 - (x - x0) ** 2) ** 0.5 for x in X2]
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

    logging.info("phi="+str(phi)+":\t"+ str(phi*180/np.pi))

    logging.info("C0="+str(C0))

    dist = np.abs(A * x0 + B * y0 + C)/(A*A+B*B)**0.5
    logging.info("\tdist="+str(dist/100000))
    logging.info("\tr="+str( r/100000))
    logging.info("\tlinear condition:")

    # r>=d
    logging.info("\t"+str((sigma3-sigma1)/2-np.abs(-1*(np.tan(phi)*(sigma1+sigma3)/2+C0)/(np.tan(phi)**2+1)**0.5) <= 0))

    q = (1 + np.sin(phi)) / (1 - np.sin(phi))

    #CoulombMohr
    # Y3 = [q*x+C0 for x in X1]
    # ax.plot(X1,Y3, marker='^')

    if (dist <= r):
        ans = "Failure occurs"
    else:
        ans = "Failure does not occur"

    logging.info("\t\t\t"+ans)

    fig.text(
        0.4, 0.1,
        ans,
        horizontalalignment='left',
        fontsize=15,
        transform=ax.transAxes
    )

    pl.savefig(output_dir+"/failure_qt/" + str(id) + ".png", dpi=200)
    # pl.savefig("./output/failure_qt/" + str(id) + ".png", dpi=200)
    # pl.show()
    pl.close()
def ro_plot_z(z_mas, ro_mas, layers_data, criteria_list_names):

    # logging.info(z_mas)
    # logging.info(ro1_vector)
    import matplotlib.ticker
    fig = pl.figure()
    # ax = fig.add_subplot(111)
    ax = fig.add_subplot(111, aspect='equal')
    pl.title(r'$\rho(z)$', size=14)
    pl.xlabel(str(RO_MIN)+"\t"+r'$\rho$'+"\t"+str(RO_MAX), size=14)
    pl.ylabel('z', size=14)

    formatter = matplotlib.ticker.NullFormatter()
    # Установка пустого форматера для оси X
    ax.xaxis.set_major_formatter(formatter)

    # отрисовываем график по каждому критерию
    ro_mas_array = np.array(ro_mas)
    for criteria_id, criteria_name in enumerate(criteria_list_names):
        if criteria_id == 0:
            # ax.plot(ro_mas_array[:, criteria_id], z_mas, label=criteria_name, marker='^')
            # logging.info("criterion_id=0")
            ax.plot(ro_mas_array[:, criteria_id], z_mas, label=criteria_name)
        else:
            ax.plot(ro_mas_array[:, criteria_id], z_mas, label=criteria_name)



    # границы слоев
    x_vector = np.linspace(RO_MIN, RO_MAX, 10)
    for layer in layers_data:
        ax.plot(x_vector, np.linspace(-layer.z, -layer.z, 10), linestyle='--', color="gray")
        # ax.plot(ro2_vector, np.linspace(-layer.z, -layer.z, 10), linestyle='--', color="gray")

    # pl.legend(loc='upper left')

    # 'best'
    # 'upper right'
    # 'upper left'
    # 'lower left'
    # 'lower right'
    # 'right'
    # 'center left'
    # 'center right'
    # 'lower center'
    # 'upper center'
    # 'center'


    # pl.savefig("./output/ro_plot.png", dpi=200)

    pl.savefig(output_dir+"/ro_plot.png", dpi=200)
    pl.legend(loc='upper left')
    pl.savefig(output_dir+"/ro_plot_legend.png", dpi=200)

    # pl.savefig("./output/ro_plot_legend.png", dpi=200)
    # pl.show()
def plot_i_ro(i_mas,ro_mas, criteria_list_names):

    fig = pl.figure()
    ax = fig.add_subplot(111)
    # ax = fig.add_subplot(111, aspect='equal')
    pl.title(r'$\rho(i)$', size=14)
    pl.xlabel(r'i', size=14)
    pl.ylabel(r'$\rho$', size=14)

    # отрисовываем график по каждому критерию
    ro_mas_array = np.array(ro_mas)
    for criteria_id, criteria_name in enumerate(criteria_list_names):
            ax.plot([i*180/np.pi for i in i_mas],ro_mas_array[:, criteria_id], label=criteria_name)


    # границы слоев
    # x_vector = np.linspace(RO_MIN, RO_MAX, 10)
    # for layer in layers_data:
    #     ax.plot(x_vector, np.linspace(-layer.z, -layer.z, 10), linestyle='--', color="gray")
    #     # ax.plot(ro2_vector, np.linspace(-layer.z, -layer.z, 10), linestyle='--', color="gray")

    # pl.legend(loc='upper left')

    # 'best'
    # 'upper right'
    # 'upper left'
    # 'lower left'
    # 'lower right'
    # 'right'
    # 'center left'
    # 'center right'
    # 'lower center'
    # 'upper center'
    # 'center'

    pl.savefig(output_dir+"/plot_i_ro.png", dpi=200)
    # pl.savefig("./output/plot_i_ro.png", dpi=200)
    pl.legend(loc='lower right')
    pl.savefig(output_dir+"/plot_i_ro_legend.png", dpi=200)
    # pl.savefig("./output/plot_i_ro_legend.png", dpi=200)
    # pl.show()


def plot_ro_k(raw_layers_data):
    fig = pl.figure()
    # ax = fig.add_subplot(111)
    ax = fig.add_subplot(111, aspect='equal')
    import matplotlib.ticker
    formatter = matplotlib.ticker.NullFormatter()
    # Установка форматера для оси X
    ax.xaxis.set_major_formatter(formatter)

    pl.title(r'$\rho_k(z)$', size=14)
    pl.xlabel(r'$\rho_k$', size=14)
    pl.ylabel('z', size=14)

    X_ro_k = [float(ro_k) * 1000 for ro_k in raw_layers_data.ro_k]
    Y_ro_k = [-1 * z for z in raw_layers_data.z]
    ax.plot(X_ro_k, Y_ro_k)

    # pl.legend(loc='upper left')
    pl.savefig(output_dir+"/ro_k_plot.png", dpi=200)
    # pl.savefig("./output/ro_k_plot.png", dpi=200)
    # pl.show()

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

def make_form_plot(theta_grid, r_grid, cur_section_cyl, j):
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
    # pl.show()
    # pl.close()
    return fig

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
    pl.savefig(output_dir+"/stress_at_point/" + str(j) + '.png', dpi=200)
    # pl.savefig('./output/' + str(j) + '.png', dpi=200)
    # pl.legend()
    # pl.show()
    pl.close()
def plot_diagram_ro_i(ro_mas, i_mas, criteria_list_names):

    title = r'$\rho(i)$'

    fig = pl.figure()
    title_size = 20

    ro_theta_vector = np.linspace(0, 2 * np.pi, 100)
    ro_r_vector = i_mas

    ro_r_grid, ro_theta_grid = np.meshgrid(ro_r_vector, ro_theta_vector)
    ro_s_r_data = np.zeros(ro_theta_grid.shape)

    for id, i in enumerate(ro_r_vector):
        ro_s_r_data[:,id] = i

    ax = pl.subplot(111, polar=True)
    title_position_x = -0.6
    title_position_y = 1

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    pl.pcolormesh(ro_theta_grid, ro_r_grid, ro_s_r_data)
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


    # prepare_one_ax(ax, pl, ro_theta_grid, ro_r_grid, s_r_data, title, title_size, Rw, title_position_x, title_position_y)

    # fig.text(
    #     0.45, 0.1,
    #   "pw=" + str(round(pw / 100000, 2)) + " atm\n" + str(z) + " m",
    #     horizontalalignment='left',
    #     fontsize=15,
    #     transform=ax.transAxes
    # )

    fig.set_label("Ro(i) diagram")
    pl.savefig(output_dir+"/stress_at_point/ro(i)" + str(id) + '.png', dpi=200)
    # pl.savefig('./output/' + str(j) + '.png', dpi=200)
    # pl.legend()
    # pl.show()
    pl.close()

def draw_points(raw_data):
    fig = pl.figure()
    ax = fig.add_subplot(111, aspect='equal')
    import matplotlib.ticker
    formatter = matplotlib.ticker.NullFormatter()
    # Установка форматера для оси X
    ax.xaxis.set_major_formatter(formatter)

    x_vector = list(raw_data.x)
    z_vector = list(-1 * raw_data.z)
    pl.xlabel(str(min(raw_data.x)) + "      " + " x " + "      " + str(max(raw_data.x)), size=14)
    pl.ylabel('z', size=14)
    ax.plot(x_vector, z_vector, marker='^', linestyle='--')
    # ax.plot(x_vector, z_vector, marker='^', linestyle='--', markevery=(0, 5))
    text1 = pl.text(0.5, 0.44, u'Trajectory')
    ax.grid(True)  # линии вспомогательной сетки
    pl.savefig(output_dir+'/points.png', dpi=300)
    # pl.savefig('./output/points.png', dpi=300)
    # pl.show()
    pl.close()

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
    ro_plot_z(zMas, ro1_vector, ro2_vector, ro3_vector)




