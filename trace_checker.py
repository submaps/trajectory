# -*- coding: utf-8 -*-
import logging
import math
import numpy as np
import matplotlib.pyplot as pl
from failure_criteria import *
from section_test import *
from failure_criteria_graph import *
from config import *
from point import *
from layer import *
from stress_cyl import StressCyl
from stress_dec import *
from section_cyl import *
from section_test import *



def rotate_general_stress_dec_mtx(a, i, sigma_H, sigma_h, sigma_v):

    # a-angle x0x' i-angle z0z'
    lxx = np.cos(a) * np.cos(i)
    lyx = -np.sin(a)
    lzx = np.cos(a) * np.sin(i)

    lxy = np.sin(a) * np.cos(i)
    lyy = np.cos(a)
    lzy = np.sin(a) * np.sin(i)

    lxz = -np.sin(i)
    lyz = 0
    lzz = np.cos(i)

    s_x = lxx ** 2 * sigma_H + lxy ** 2 * sigma_h + lxz ** 2 * sigma_v
    s_y = lyx ** 2 * sigma_H + lyy ** 2 * sigma_h + lyz ** 2 * sigma_v
    s_z = lzx ** 2 * sigma_H + lzy ** 2 * sigma_h + lzz ** 2 * sigma_v

    t_x_y = lxx * lyx * sigma_H + lxy * lyy * sigma_h + lxz * lyz * sigma_v
    t_y_z = lyx * lzx * sigma_H + lyy * lzy * sigma_h + lyz * lzz * sigma_v
    t_x_z = lzx * lxx * sigma_H + lzy * lxy * sigma_h + lzz * lxz * sigma_v

    return StressDec(
        s_x,
        s_y,
        s_z,
        t_y_z,
        t_x_z,
        t_x_y
    )


def calculate_stress_cyl(Rw, v_fr, pw, r, theta, s_0):
    # s_x_0, s_y_0, s_z_0, t_yz_0, t_xz_0, t_xy_0
    s_x_0 = s_0.s_x
    s_y_0 = s_0.s_y
    s_z_0 = s_0.s_z
    t_yz_0 = s_0.t_y_z
    t_xz_0 = s_0.t_x_z
    t_xy_0 = s_0.t_x_y

    s_r = (s_x_0 + s_y_0) / 2 * (1 - Rw ** 2 / r ** 2) + (s_x_0 - s_y_0) / 2 * (
           1 + 3 * Rw ** 4 / r ** 4 - 4 * Rw ** 2 / r ** 2) * np.cos(2 * theta) + \
          +t_xy_0 * (
           1 + 3 * Rw ** 4 / r ** 4 - 4 * Rw ** 2 / r ** 2) * np.sin(2 * theta) \
          - pw * Rw ** 2 / r ** 2
    s_theta = (s_x_0 + s_y_0) / 2 * (1 + Rw ** 2 / r ** 2) - (s_x_0 - s_y_0) / 2 * (
        1 + 3 * Rw ** 4 / r ** 4) * np.cos(2 * theta)
    - t_xy_0 * (1 + 3 * Rw ** 4 / r ** 4) * np.sin(2 * theta) + pw * Rw ** 2 / r ** 2
    # s_x_0+s_y_0-2*(s_x_0-s_y_0)*np.cos(2*theta)-t_xy_0*np.sin(2*theta)-pw
    # (s_x_0+s_y_0)/2*(1+Rw**2/r**2)-(s_x_0-s_y_0)/2*(1+3*Rw**4/r**4)*np.cos(2*theta)-t_xy_0*(1+3*Rw**4/r**4)*sin(2*theta)-pw*Rw**2/r**2

    s_z = s_z_0 - v_fr * (
        2 * (s_x_0 - s_y_0) * Rw ** 2 / r ** 2 * np.cos(2 * theta) +
        4 * t_xy_0 * Rw ** 2 / r ** 2 * np.sin(2 * theta))

    # s_z = s_z_0 - v_fr * (
    #     2 * (s_x_0 - s_y_0) * Rw ** 2 / r ** 2 * np.cos(2 * theta) +
    #     4 * t_xy_0 * Rw ** 2 / r ** 2 * np.sin(2 * theta))

    #! s_y_0-s_x_0
    t_r_theta = (s_x_0 - s_y_0) / 2 * (1 - 3 * Rw ** 4 / r ** 4 + 2 * Rw ** 2 / r ** 2) * np.sin(
        2 * theta) + t_xy_0 * (1 - 3 * Rw ** 4 / r ** 4 + 2 * Rw ** 2 / r ** 2) * np.cos(2 * theta)

    t_theta_z = (-1*t_xz_0 * np.sin(theta) + t_yz_0 * np.cos(theta)) * (1 + Rw ** 2 / r ** 2)
    t_r_z = (t_xz_0 * np.cos(theta) + t_yz_0 * np.sin(theta)) * (1 - Rw ** 2 / r ** 2)

    return StressCyl(s_r,
                      s_theta,
                      s_z,
                      t_r_theta,
                      t_theta_z,
                      t_r_z)


def calculate_stress_cyl_at_well(nu, pw, theta, s_0):
    # s_x_0, s_y_0, s_z_0, t_yz_0, t_xz_0, t_xy_0
    s_x_0 = float(s_0.s_x)
    s_y_0 = float(s_0.s_y)
    s_z_0 = float(s_0.s_z)
    t_yz_0 = float(s_0.t_y_z)
    t_xz_0 = float(s_0.t_x_z)
    t_xy_0 = float(s_0.t_x_y)

    # s_r = float(pw)
    # s_theta = s_x_0 + s_y_0 - 2 * (s_x_0 - s_y_0) * np.cos(2 * theta) - t_xy_0 * np.sin(2 * theta) + pw
    # s_z = s_z_0 - nu * (2 * (s_x_0 - s_y_0) * np.cos(2 * theta) + 4 * t_xy_0 * np.sin(2 * theta))
    # t_r_theta = 0
    # t_theta_z = 2 * (-1*t_xz_0 * np.sin(theta) + t_yz_0 * np.cos(theta))
    # t_r_z = 0


    s_r = -pw
    s_theta = s_x_0 + s_y_0 - 2 * (s_x_0 - s_y_0) * np.cos(2 * theta) - 4 * t_xy_0 * np.sin(2 * theta) + pw
    s_z = s_z_0 - nu * (2 * (s_x_0 - s_y_0) * np.cos(2 * theta) + 4 * t_xy_0 * np.sin(2 * theta))
    t_r_theta = 0
    t_theta_z = 2 * (-t_xz_0 * np.sin(theta) + t_yz_0 * np.cos(theta))
    t_r_z = 0

    # s_r, s_theta, s_z, t_r_theta, t_theta_z, t_r_z

    stressCyl = StressCyl(
                    s_r,
                    s_theta,
                    s_z,
                    t_r_theta,
                    t_theta_z,
                    t_r_z)
    return stressCyl


def load_points(source_file):
    # import pandas
    # data = pandas.read_csv(source_file, ';', decimal=',')
    raw_data = load_trajectory(source_file)
    # id;x;y;z
    pointsData = []
    for (id, x, y, z) in zip(raw_data.id, raw_data.x, raw_data.y, raw_data.z):
        t = Point(id, x, y, z)
        t.printPoint()
        pointsData.append(t)
    logging.info("________________________")
    return pointsData, raw_data


def load_trajectory(source_file):
    import pandas
    data = pandas.read_csv(source_file, ';', decimal=',')
    # id;x;y;z
    return data


def load_layers(source_file):
    import pandas
    raw_layers_data = pandas.read_csv(source_file, ';', decimal=',')
    # id;z;h_k;ro_k;nu_k
    # data.id,data.z,data.h_k,data.ro_k,data.nu_k
    layersData = []
    for (id, z, h_k, ro_k, nu_k, phi_k, C0_k) in zip(raw_layers_data.id, raw_layers_data.z, raw_layers_data.h_k, raw_layers_data.ro_k, raw_layers_data.nu_k, raw_layers_data.phi_k,
                                                     raw_layers_data.C0_k):
        l = Layer(id, z, h_k, ro_k, nu_k, float(phi_k) * np.pi / 180, float(C0_k) * 1000000)
        l.printLayer()
        layersData.append(l)
    logging.info("________________________")
    return layersData, raw_layers_data


def get_layer_from_z(z, layersData):
    for l in layersData:
        if z <= l.z:
            return l
    logging.error("Layer not found!" + str(z))
    # print "Layer not found!" + str(z)
    return layersData[0]


def calculate_section(j, i, a, pw, sigma_H, sigma_h, sigma_v, r_vector, theta_vector, gridShape, Rw, z, nu):
    # главное напряжение в данной точке
    s_0 = rotate_general_stress_dec_mtx(a, i, sigma_H, sigma_h, sigma_v)
    logging.info("printStressDec")
    logging.info("trace comparison: "+str(s_0.getTrace()-(sigma_H+sigma_h+sigma_v)))
    s_0.printStressDec()

    # массивы для каждого компонента напряжения для каждой точки
    # размерности gridShape как r_grid, theta_grid
    s_r_data = np.zeros(gridShape)
    s_theta_data = np.zeros(gridShape)
    s_z_data = np.zeros(gridShape)

    t_r_theta_data = np.zeros(gridShape)
    t_theta_z_data = np.zeros(gridShape)
    t_r_z_data = np.zeros(gridShape)

    currentGeneralStressListInPoint = []
    # массив StressCyl объектов для проверки критериев
    stressAtWell = []

    # проходим сначала по радиусам, затем по углам
    # rc индекс текущего радиуса, r текущий радиус
    # thetac индекс текущего угла, theta текущий угол


    for rc, r in enumerate(r_vector):
        for thetac, theta in enumerate(theta_vector):
            # значения тензора напряжения в этой точке
            currentStressCyl = calculate_stress_cyl(Rw, nu, pw, r, theta, s_0)
            if r == Rw:
                stressAtWell.append(currentStressCyl)
            s_r_data[thetac][rc] = currentStressCyl.s_r
            s_theta_data[thetac][rc] = currentStressCyl.s_theta
            s_z_data[thetac][rc] = currentStressCyl.s_z
            t_r_theta_data[thetac][rc] = currentStressCyl.t_r_theta
            t_theta_z_data[thetac][rc] = currentStressCyl.t_theta_z
            t_r_z_data[thetac][rc] = currentStressCyl.t_r_z

    # тип StressCylMas в котором находятся массивы s_r,s_theta,...
    # для текущего сечения
    curSectionCyl = SectionCyl(s_r_data,
                               s_theta_data,
                               s_z_data,
                               t_r_theta_data,
                               t_theta_z_data,
                               t_r_z_data,
                               i, a, pw, Rw, z, stressAtWell)
    return curSectionCyl


def calculate_section_at_well(j, i, a, pw, sigma_H, sigma_h, sigma_v, r_vector, theta_vector, gridShape, Rw, z, nu):
    # главное напряжение в данной точке
    s_0 = rotate_general_stress_dec_mtx(a, i, sigma_H, sigma_h, sigma_v)
    # logging.info("stress trace comparison: "+str(s_0.getTrace()-(sigma_H+sigma_h+sigma_v)))
    # s_0.printStressDec()
    stress_at_well = []
    for theta in theta_vector:
        currentStressCyl = calculate_stress_cyl_at_well(nu, pw, theta, s_0)
        stress_at_well.append(currentStressCyl)

    return stress_at_well

def initialize_logger():

    logger = logging.getLogger()
    logger.setLevel(log_level)
    # create console handler and set level to info
    handler = logging.StreamHandler()
    handler.setLevel(log_level)
    formatter = logging.Formatter("%(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # create error file handler and set level to error
    handler = logging.FileHandler("experiment/trace_checker_error.log")
    handler.setLevel(logging.ERROR)
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # create debug file handler and set level to info
    handler = logging.FileHandler("experiment/trace_checker.log")
    handler.setLevel(log_level)
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def calculate_angle(point_prev, point):
    # i # угол z0z' в радианах
    # a # угол x0x' в радианах
    # находим углы наклона по координатам текущей и предыдущей точки

    i = np.pi/2-math.atan2(point.z-point_prev.z,point.x-point_prev.x)

    if point == point_prev:
        i = 0
        a = 0

    a = 0


    return i, a


def print_inf(point, point_prev, l, a, i, sigma_H, sigma_h, sigma_v):
    logging.info(str(int(point.id)) + ":\t\tz: " + str(point.z) + " l:" + str(l.id) + ";\t i_deg: " \
                 + str(round(i * 180 / np.pi, 3)) + "; a_deg: " + str(round(a * 180 / np.pi, 3)))
    logging.info("l.z:" + str(l.z) + " h:" + str(l.h_k) + " z:" + str(point.z) + " dz=" + str(point.z - point_prev.z))
    logging.info("sigma_H: " + str(sigma_H / 100000) + " atm \t sigma_h: " + str(
        sigma_h / 100000) + " atm \t sigma_v: " + str(sigma_v / 100000) \
                 + "\tpw: undef" + "\tro: undef")
    logging.info("s_0 after rotate:")
    rotatedStressMtx = rotate_general_stress_dec_mtx(a, i, sigma_H, sigma_h, sigma_v)
    rotatedStressMtx.printStressDec()
    # logging.info("s_general:")
    # logging.info(rotatedStressMtx.getGeneralStress())


def calculate_trajectory():
    # задаем настройки для логгера
    initialize_logger()
    logging.info("________________________\n" + "calculate each section" + "\n________________________")
    # pointsData массив объектов типа Point
    # raw_data набор значений из csv
    points_data, raw_points_data = load_points(input_trajectory)
    layers_data, raw_layers_data = load_layers(input_layers)

    # отрисовка траектории
    draw_points(raw_points_data)

    # сетка по r и theta
    r_grid, theta_grid = np.meshgrid(r_vector, theta_vector)

    # размерность сетки # 100, 25
    gridShape = r_grid.shape
    logging.info("r_grid.shape=" + str(r_grid.shape))

    sections_mas = []
    # generalMasPoint = []
    z_mas = []  # массив координат точек по z
    i_mas = []  # массив углов

    from config import ro_init
    sigma_v = 0
    sigma_H = 0

    point_prev = points_data[0]  # в качестве первой prev точки берем первую из массива
    # проходим по всем точкам, загруженным из csv
    for point in points_data:
        logging.info("----------------------------------------------------")
        l = get_layer_from_z(point.z, layers_data)  # находим какому слою принадлежит эта точка
        # sigma_H, sigma_h, sigma_v главные напряжения в породе по z
        # l.ro_k плотность породы в слое, в файле в г/см3 переводим в кг/м3
        pw = ro_init*g*point.z
        sigma_v += - g * float(l.ro_k * 1000) * (point.z - point_prev.z)

        # sigma_v = l.ro_k*g*point.z
        sigma_H = l.nu_k / (1.0 - l.nu_k) * sigma_v
        sigma_h = sigma_H  # считаем, что напряжения распора равны

        # a # угол x0x' в радианах
        # i # угол z0z' в радианах
        # находим углы наклона по координатам текущей и предыдущей точки
        i, a = calculate_angle(point_prev, point)
        print_inf(point, point_prev, l, a, i, sigma_H, sigma_h, sigma_v)

        # расчет поля напряжений в текущем сечении
        cur_section_cyl = calculate_section(point.id, i, a, pw, sigma_H, sigma_h, sigma_v, r_vector, theta_vector, gridShape, Rw,
                                         point.z, l.nu_k)
        sections_mas.append(cur_section_cyl)
        # добавляем данные (сечения) для построения графиков в массив

        i_mas.append(i)  # массив углов i
        # pwMas.append(pw)  # массив давлений по высоте
        z_mas.append(-1 * point.z)  # массив по z
        point_prev = point

    logging.info(z_mas)
    # logging.info(ro1_vector)

    for j in range(0, sections_mas.__len__()):
        print "img ", j
        make_stress_plot(
            theta_grid,
            r_grid,
            sections_mas[j],
            j
        )


def ro_analysis():

    # задаем настройки для логгера
    initialize_logger()
    logging.info("________________________\n" + "Trace checker ro_analysis starting" + "\n________________________")
    # pointsData массив объектов типа Point
    # raw_data набор значений из csv
    points_data, raw_points_data = load_points(input_trajectory)
    layers_data, raw_layers_data = load_layers(input_layers)

    # отрисовка траектории
    draw_points(raw_points_data)

    # сетка по r и theta
    r_grid, theta_grid = np.meshgrid(r_vector, theta_vector)

    # размерность сетки # 100, 25
    gridShape = r_grid.shape
    logging.info("r_grid.shape=" + str(r_grid.shape))

    z_mas = []  # массив координат точек по z
    i_mas = []  # массив углов
    ro_mas = []  # минимальные плотности по каждому критерии

    # pwMas = []  # массив давления
    # printMas = []  # массив для печати
    point_prev = points_data[0]  # в качестве первой prev точки берем первую из массива
    sigma_v = 0
    # проходим по всем точкам, загруженным из csv
    for point in points_data:
        logging.info("----------------------------------------------------")
        l = get_layer_from_z(point.z, layers_data)  # находим какому слою принадлежит эта точка
        # sigma_H, sigma_h, sigma_v главные напряжения в породе по z
        # l.ro_k плотность породы в слое, в файле в г/см3 переводим в кг/м3

        sigma_v += - g * float(l.ro_k * 1000) * (point.z - point_prev.z)

        # sigma_v = l.ro_k*g*point.z
        sigma_H = l.nu_k / (1.0 - l.nu_k) * sigma_v
        sigma_h = sigma_H  # считаем, что напряжения распора равны

        # a # угол x0x' в радианах
        # i # угол z0z' в радианах
        # находим углы наклона по координатам текущей и предыдущей точек
        i, a = calculate_angle(point_prev, point)

        print_inf(point, point_prev, l, a, i, sigma_H, sigma_h, sigma_v)

        # минимальная плотность раствора в текущем сечении
        # возвращается list из ro каждого критерия
        ro_mas.append(get_min_ro_in_section(point.id, i, a, sigma_H, sigma_h, sigma_v,
                                            r_vector, theta_vector, gridShape,
                                            Rw, point.z, l.nu_k, l.phi_k, l.C0_k))

        #[800 800 800 800 800 1300 1300 1300 1300]

        i_mas.append(i)  # массив углов i
        z_mas.append(-1*point.z)  # массив по z
        point_prev = point

    criteria_list_names = FailureCriteria(RO_MAX, "test").criteria_names

    for z, i, ro in zip(z_mas, i_mas, ro_mas):
        ro_str=""
        for ro_id in ro:
            ro_str += str(ro_id)+";"
        logging.info(str(z)+";\ti: ;"+str(i)+"\t"+";\ti_dec;"+str(round(i*180/np.pi,2))+";\t"+ro_str)

    ro_plot_z(z_mas,
              ro_mas,
              layers_data,
              criteria_list_names
              )

    plot_ro_k(raw_layers_data)

    plot_i_ro(i_mas,ro_mas, criteria_list_names)

    plot_diagram_ro_i(ro_mas, i_mas, criteria_list_names)

if __name__ == '__main__':
    if target == "ro_analysis":
        ro_analysis()

    if target == "calculate_trajectory":
        calculate_trajectory()


