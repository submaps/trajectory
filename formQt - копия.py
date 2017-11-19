# -*- coding:utf8 -*-
import sys
import numpy as np
from PyQt4 import QtGui
from PyQt4 import QtCore

import logging

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg \
    import FigureCanvasQTAgg as FigureCanvas


from failure_criteria_graph import *
from config import *
from section_test import *
from trace_checker import calculate_section_at_well
from stress_cyl import *

class Monitor(FigureCanvas, QtCore.QObject):
    def __init__(self, cur_section_cyl):
        self.fig = Figure()
        # initialize the figure canvas
        FigureCanvas.__init__(self, self.fig)
        self.drawPlot(cur_section_cyl)

        # self.fig = Figure(facecolor='#E8D6BB', edgecolor='#000000')
        # self.fig = Figure()
        # self.ax = self.fig.add_subplot(111)
        #
        # # initialize the figure canvas
        # FigureCanvas.__init__(self, self.fig)
        #
        # # set up the display limits for the figure
        # self.ax.set_xlim(0, 30)
        # self.ax.set_ylim(0, 1)
        #
        # # turn off autoscaling
        # self.ax.set_autoscale_on(False)
        # line, = self.ax.plot(np.random.rand(100), 'o', picker=5)
        # self.fig.canvas.draw()
        # self.stopButtonAction = False

    def drawPlot(self, cur_section_cyl):
        self.cur_section_cyl = cur_section_cyl


        r_grid, theta_grid = np.meshgrid(r_vector, theta_vector)
        # размерность сетки # 100, 25

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

        title_size = 20

        title_position_x = -0.6
        title_position_y = 1

        ax1 = self.fig.add_subplot(331, polar=True)
        self.prepare_one_ax_for_canvas(ax1, theta_grid, r_grid, s_r_data, title1, title_size, Rw,
                                       title_position_x,
                                       title_position_y)

        ax2 = self.fig.add_subplot(334, polar=True)
        self.prepare_one_ax_for_canvas(ax2, theta_grid, r_grid, s_theta_data, title2, title_size, Rw,
                                       title_position_x,
                                       title_position_y)

        ax3 = self.fig.add_subplot(337, polar=True)
        self.prepare_one_ax_for_canvas(ax3, theta_grid, r_grid, s_z_data, title3, title_size, Rw,
                                       title_position_x,
                                       title_position_y)

        ax4 = self.fig.add_subplot(333, polar=True)
        self.prepare_one_ax_for_canvas(ax4,  theta_grid, r_grid, t_r_theta_data, title4, title_size, Rw,
                                       title_position_x,
                                       title_position_y)

        ax5 = self.fig.add_subplot(336, polar=True)
        self.prepare_one_ax_for_canvas(ax5, theta_grid, r_grid, t_theta_z_data, title5, title_size, Rw,
                                       title_position_x,
                                       title_position_y)

        ax6 = self.fig.add_subplot(339, polar=True)
        self.prepare_one_ax_for_canvas(ax6, theta_grid, r_grid, t_r_z_data, title6, title_size, Rw,
                                       title_position_x,
                                       title_position_y)

        self.fig.text(
            0.45, 0.1,
            "i=" + str(round(i * 180 / np.pi, 2)) + "\npw=" + str(round(pw / 100000, 2)) + " atm\n" + str(z) + " m",
            horizontalalignment='left',
            fontsize=15,
            transform=ax1.transAxes
        )

        self.fig.canvas.draw()
        self.fig.set_label("Main graph")

    def prepare_one_ax_for_canvas(self, ax, theta_grid, r_grid, data, title, title_size, Rw, title_position_x,
                                  title_position_y):
        # color_mesh = fig.matrix_axes.pcolormesh(data, cmap=color_map)
        ##
        # add tick labels
        # self.matrix_axes.set_yticklabels(labels_y)
        # self.matrix_axes.set_xticklabels(labels_x)
        # self.figure.autofmt_xdate(rotation=30)  # rotate x axis labels to fit more
        # plot color bar
        # colorbar = self.figure.colorbar(color_mesh, cax=self.color_axes, orientation='vertical')
        ##
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        # color_mesh = self.matrix_axes.pcolormesh(data, cmap=color_map)
        color_mesh = ax.pcolormesh(theta_grid, r_grid, data)
        # self.pcolormesh(theta_grid, r_grid, data)
        colorbar = self.figure.colorbar(color_mesh)

        # colorbar = self.figure.colorbar(color_mesh, cax=self.color_axes, orientation='vertical')
        # self.colorbar(pad=0.1)
        # pl.clim(minval, maxval)
        ax.set_thetagrids(np.array([0, 90, 180, 270]), ['0', '90', '180', '270'])
        # ax.set_thetagrids(np.array([0, 90, 180, 270]), ['0', '90', '180', '270'],fontsize=8)
        # ax.set_rgrids(radii=[Rw * 1, Rw * 2], labels=['1', ' 2 '], angle=90, fontsize=15)
        # ax.set_rgrids(radii=[Rw * 1, Rw * 2], labels=['$R_w$', '$2R_w$'], angle=90, fontsize=8)
        ax.set_rgrids(radii=[Rw * 1, Rw * 2], angle=90, fontsize=8)
        # ax.set_rgrids(radii=[Rw * 1, Rw * 2], labels=['', ''], angle=90, fontsize=8)

        ax.grid(True, color='black', linestyle='-', linewidth=1, axis='y')
        self.fig.text(title_position_x, title_position_y,
                title,
                horizontalalignment='left',
                fontsize = title_size,
                transform = ax.transAxes)
        ax.set_title(title,loc='left')

    def updatePlot(self, cur_section_cyl):
        self.fig.clear()
        self.drawPlot(cur_section_cyl)
        self.fig.canvas.draw()

    @QtCore.pyqtSignature("")
    def zoomIn(self):
        """
        Увеличиваем в 2 раза
        """
        self.stopButtonAction = True
        while self.stopButtonAction:

            start, end = self.ax.get_xaxis().get_view_interval()
            print start, end, end - start
            if (end - start > 0.0001):
                self.ax.get_xaxis().set_view_interval(start + ((end - start) / 4)
                                                      , end - ((end - start) / 4), True)

            print "Before draw"
            self.fig.canvas.draw()
            print "Before Qt"
            QtGui.qApp.processEvents();

    @QtCore.pyqtSignature("")
    def zoomOut(self):
        """
        Увеньшает в 2 раза
        """
        self.stopButtonAction = True
        while self.stopButtonAction:
            start, end = self.ax.get_xaxis().get_view_interval()
            self.ax.get_xaxis().set_view_interval(start - ((end - start) / 2)
                                                  , end + ((end - start) / 2), True)

            self.fig.canvas.draw()
            QtGui.qApp.processEvents();

    @QtCore.pyqtSignature("")
    def panLeft(self):
        """
        Сдвигаем на 1/10
        """
        self.stopButtonAction = True
        while self.stopButtonAction:
            start, end = self.ax.get_xaxis().get_view_interval()
            interval = (end - start) / 10
            self.ax.get_xaxis().set_view_interval(start - interval
                                                  , end - interval, True)

            self.fig.canvas.draw()
            QtGui.qApp.processEvents();

    @QtCore.pyqtSignature("")
    def panRight(self):
        self.stopButtonAction = True
        while self.stopButtonAction:
            start, end = self.ax.get_xaxis().get_view_interval()
            interval = (end - start) / 10
            self.ax.get_xaxis().set_view_interval(start + interval
                                                  , end + interval, True)
            self.fig.canvas.draw()
            QtGui.qApp.processEvents();

    @QtCore.pyqtSignature("")
    def stopAction(self):
        self.stopButtonAction = False


class MyWindow(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        # initial consts
        self.window_initialize_logger()
        self.Rw = 0.01
        self.angle_i_dec = 45
        self.z = 1000
        self.ro = 1300
        self.ro_k = 2300
        self.nu_k = 0.35
        self.phi_k_dec = 44
        self.C0_k = 25000000
        self.pw = self.ro * 9.81 * self.z
        self.sigma_v = stress_sign*self.ro_k * 9.81 * self.z
        self.sigma_H = self.nu_k / (1 - self.nu_k) * self.sigma_v
        self.fcriteria = FailureCriteria(RO_MIN, 777)
        r_grid, theta_grid = np.meshgrid(r_vector, theta_vector)
        # размерность сетки # 100, 25
        gridShape = r_grid.shape

        self.cur_section_cyl = self.window_calculate_section(id, (self.angle_i_dec * np.pi / 180), 0, self.pw,
                                            self.sigma_H, self.sigma_H, self.sigma_v,
                                            r_vector, theta_vector, gridShape, Rw,
                                            self.z, self.nu_k)
        # self.cur_section_cyl = cur_section_cyl
        self.m = Monitor(self.cur_section_cyl)
        self.create_ui()
        self.updateView()

    @QtCore.pyqtSignature("")
    def startAction(self):
        logging.info("============Start action============")
        self.updateView()
        logging.info("z: "+str(self.z))
        logging.info("pw: "+str(self.pw))
        logging.info("i : " + str(self.angle_i_dec))
        logging.info("ro: "+str(self.ro))

        id = 777
        r_grid, theta_grid = np.meshgrid(r_vector, theta_vector)
        # размерность сетки # 100, 25
        gridShape = r_grid.shape

        self.cur_section_cyl = self.window_calculate_section(id, self.angle_i_dec * np.pi / 180, 0, self.pw,
                                            self.sigma_H, self.sigma_H, self.sigma_v,
                                            r_vector, theta_vector, gridShape, Rw,
                                            self.z, self.nu_k)

        self.m.updatePlot(self.cur_section_cyl)

    def window_initialize_logger(self):
        logger = logging.getLogger()
        logger.setLevel(log_level)
        # create console handler and set level to info
        handler = logging.StreamHandler()
        handler.setLevel(log_level)
        formatter = logging.Formatter("%(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)

        # create error file handler and set level to error
        handler = logging.FileHandler("output/form_qt_error.log")
        handler.setLevel(logging.ERROR)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

        # create debug file handler and set level to info
        handler = logging.FileHandler("output/form_qt.log")
        handler.setLevel(log_level)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    def window_calculate_section(self, j, i, a, pw, sigma_H, sigma_h, sigma_v, r_vector, theta_vector, gridShape, Rw, z, nu):
        # главное напряжение в данной точке
        logging.info("s_0:")
        s_0 = self.window_rotate_general_stress_dec_mtx(a, i, sigma_H, sigma_h, sigma_v)
        logging.info(s_0.getMtx33())
        logging.info("\ttrace comparison: "+str(s_0.getTrace()-(sigma_H+sigma_h+sigma_v)))
        # logging.info("printStressDec")

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

    def window_rotate_general_stress_dec_mtx(self, a, i, sigma_H, sigma_h, sigma_v):

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

        s_x     = lxx ** 2 * sigma_H + lxy ** 2 * sigma_h + lxz ** 2 * sigma_v
        s_y     = lyx ** 2 * sigma_H + lyy ** 2 * sigma_h + lyz ** 2 * sigma_v
        s_z     = lzx ** 2 * sigma_H + lzy ** 2 * sigma_h + lzz ** 2 * sigma_v

        t_x_y = lxx * lyx * sigma_H + lxy * lyy * sigma_h + lxz * lyz * sigma_v
        t_y_z = lyx * lzx * sigma_H + lyy * lzy * sigma_h + lyz * lzz * sigma_v
        t_x_z = lzx * lxx * sigma_H + lzy * lxy * sigma_h + lzz * lxz * sigma_v

        #a,b
        # a = a
        # b = i
        # rotate_mtx = np.matrix([
        #             [np.sin(b)**2,np.cos(b)**2 * np.cos(a)**2,np.cos(b)**2 * np.sin(a)**2],
        #             [0,np.sin(a)**2,np.cos(a)**2],
        #             [np.cos(b)**2,np.sin(b)**2 * np.cos(a)**2,np.sin(b)**2*np.sin(a)**2],
        #             [0, -np.sin(a)*np.cos(a)* np.sin(b), np.sin(a)*np.cos(a)*np.sin(b)],
        #             [-np.sin(b)*np.cos(b),np.sin(b)*np.cos(b)*np.cos(a)**2,np.sin(b)*np.cos(b)*np.sin(a)**2],
        #             [0,-np.sin(a)*np.cos(a)*np.cos(b),np.sin(a)*np.cos(a)*np.cos(b)]
        # ])
        #
        # [[s_x], [s_y], [s_z], [t_x_y], [t_y_z], [t_x_z]] = rotate_mtx*[[sigma_v],[sigma_H],[sigma_h]]

        return StressDec(
            s_x,
            s_y,
            s_z,
            t_y_z,
            t_x_z,
            t_x_y
        )

    def window_calculate_stress_cyl_at_well(self,nu, pw, theta, s_0):
        # s_x_0, s_y_0, s_z_0, t_yz_0, t_xz_0, t_xy_0
        s_x_0 = float(s_0.s_x)
        s_y_0 = float(s_0.s_y)
        s_z_0 = float(s_0.s_z)
        t_yz_0 = float(s_0.t_y_z)
        t_xz_0 = float(s_0.t_x_z)
        t_xy_0 = float(s_0.t_x_y)

        s_r = float(pw)
        s_theta = s_x_0 + s_y_0 - 2 * (s_x_0 - s_y_0) * np.cos(2 * theta) - t_xy_0 * np.sin(2 * theta) - pw
        s_z = s_z_0 - nu * (2 * (s_x_0 - s_y_0) * np.cos(2 * theta) + 4 * t_xy_0 * np.sin(2 * theta))
        t_r_theta = 0
        t_theta_z = 2 * (-1 * t_xz_0 * np.sin(theta) + t_yz_0 * np.cos(theta))
        t_r_z = 0

        # s_r, s_theta, s_z, t_r_theta, t_theta_z, t_r_z

        stressCyl = stress_cyl(stress_sign * s_r,
                               stress_sign * s_theta,
                               stress_sign * s_z,
                               stress_sign * t_r_theta,
                               stress_sign * t_theta_z,
                               stress_sign * t_r_z)
        return stressCyl

    def window_calculate_section_at_well(self, j, i, a, pw, sigma_H, sigma_h, sigma_v, r_vector, theta_vector, gridShape, Rw, z, nu):
        # главное напряжение в данной точке
        s_0 = self.window_rotate_general_stress_dec_mtx(a, i, sigma_H, sigma_h, sigma_v)
        # s_0.printStressDec()

        # массивы для каждого компонента напряжения для каждой точки
        # размерности gridShape как r_grid, theta_grid
        # s_r_data = np.zeros(gridShape)
        # s_theta_data = np.zeros(gridShape)
        # s_z_data = np.zeros(gridShape)
        #
        # t_r_theta_data = np.zeros(gridShape)
        # t_theta_z_data = np.zeros(gridShape)
        # t_r_z_data = np.zeros(gridShape)

        currentGeneralStressListInPoint = []
        # массив StressCyl объектов для проверки критериев
        stress_at_well = []

        for theta in theta_vector:
            currentStressCyl = self.window_calculate_stress_cyl_at_well(nu, pw, theta, s_0)
            stress_at_well.append(currentStressCyl)

        # тип StressCylMas в котором находятся массивы s_r,s_theta,...
        # для текущего сечения
        # curSectionCyl = SectionCyl(s_r_data,
        #                            s_theta_data,
        #                            s_z_data,
        #                            t_r_theta_data,
        #                            t_theta_z_data,
        #                            t_r_z_data,
        #                            i, a, pw, Rw, z, stress_at_well)
        # return curSectionCyl
        return stress_at_well

    def updateView(self):
        self.angle_i_dec = float(self.angle_i_dec_text.toPlainText())
        self.ro_k = float(self.ro_k_text.toPlainText())
        self.ro = float(self.ro_text.toPlainText())
        self.z = float(self.z_text.toPlainText())

        if self.ro_k_text.toPlainText() == "":
            self.sigma_v = float(self.sigma_v_text.toPlainText())
            self.sigma_H = float(self.sigma_H_text.toPlainText())
        else:
            # TO DO stress_sign?
            self.sigma_v = stress_sign*self.ro_k * 9.81 * self.z
            self.sigma_H = self.nu_k / (1 - self.nu_k) * self.sigma_v
        # плотность бурового раствора
        if not self.ro_text.toPlainText() == "":
            self.pw = self.ro * 9.81 * self.z
        else:
            self.pw = float(self.pw_text.toPlainText())

        self.phi_k_dec = float(self.phi_k_text.toPlainText())
        self.C0_k = float(self.C0_k_text.toPlainText())

        self.pw_text.setPlainText(str(self.pw))
        self.sigma_v_text.setPlainText(str(self.sigma_v))
        self.sigma_H_text.setPlainText(str(self.sigma_H))

    def create_ui(self):
        self.setWindowTitle(u'Stress plots in section')
        # app = QtGui.QApplication(sys.argv)

        self.startButton = QtGui.QPushButton("Start")
        self.checkButton = QtGui.QPushButton("check failure")
        # self.checkButton.clicked.connect(self.checkAction)


        self.angle_i_dec_label = QtGui.QLabel("i = ")
        self.angle_i_dec_text = QtGui.QPlainTextEdit()
        self.angle_i_dec_text.setPlainText(str(self.angle_i_dec))

        self.angle_i_dec_text.setFixedHeight(25)
        self.angle_i_dec_text.setFixedWidth(100)

        self.sigma_v_label = QtGui.QLabel("\tsigma_v = ")
        self.sigma_H_label = QtGui.QLabel("\tsigma_H = ")

        self.sigma_v_text = QtGui.QPlainTextEdit()
        self.sigma_v_text.setPlainText(str(self.sigma_v))
        self.sigma_v_text.setFixedHeight(25)
        self.sigma_v_text.setFixedWidth(150)

        self.sigma_H_text = QtGui.QPlainTextEdit()
        self.sigma_H_text.setPlainText(str(self.sigma_H))
        self.sigma_H_text.setFixedHeight(25)
        self.sigma_H_text.setFixedWidth(150)

        self.ro_k_label = QtGui.QLabel("\tro_k = ")

        self.ro_k_text = QtGui.QPlainTextEdit()
        self.ro_k_text.setPlainText(str(self.ro_k))
        self.ro_k_text.setFixedHeight(25)
        self.ro_k_text.setFixedWidth(150)

        self.phi_k_label = QtGui.QLabel("\tphi_k = ")

        self.phi_k_text= QtGui.QPlainTextEdit()
        self.phi_k_text.setPlainText(str(self.phi_k_dec))
        self.phi_k_text.setFixedHeight(25)
        self.phi_k_text.setFixedWidth(150)

        self.C0_k_label = QtGui.QLabel("\tC0_k = ")

        self.C0_k_text= QtGui.QPlainTextEdit()
        self.C0_k_text.setPlainText(str(self.C0_k))
        self.C0_k_text.setFixedHeight(25)
        self.C0_k_text.setFixedWidth(150)

        self.pw_label = QtGui.QLabel("\tpw = ")

        self.pw_text = QtGui.QPlainTextEdit()
        self.pw_text.setPlainText(str(self.pw))
        self.pw_text.setFixedHeight(25)
        self.pw_text.setFixedWidth(150)

        self.z_label = QtGui.QLabel("z = ")
        self.z_text = QtGui.QPlainTextEdit()
        self.z_text.setFixedHeight(25)
        self.z_text.setFixedWidth(35)
        self.z_text.setPlainText(str(self.z))

        self.ro_label = QtGui.QLabel("ro = ")
        self.ro_text = QtGui.QPlainTextEdit()
        self.ro_text.setFixedHeight(25)
        self.ro_text.setFixedWidth(35)
        self.ro_text.setPlainText(str(self.ro))

        self.scrollbar = QtGui.QScrollBar()

        # layout = QtGui.QVBoxLayout()
        layout = QtGui.QGridLayout()
        layout.setSpacing(10)
        layout.addWidget(self.startButton, 0, 0)
        layout.addWidget(self.checkButton, 0, 1)

        layout.addWidget(self.angle_i_dec_label, 1, 0)
        layout.addWidget(self.angle_i_dec_text, 1, 1)

        layout.addWidget(self.sigma_v_label, 1, 2)
        layout.addWidget(self.sigma_H_label, 2, 2)

        layout.addWidget(self.sigma_v_text, 1, 3)
        layout.addWidget(self.sigma_H_text, 2, 3)

        layout.addWidget(self.ro_k_label, 1, 4)
        layout.addWidget(self.ro_k_text, 1, 5)

        layout.addWidget(self.phi_k_label, 2, 4)
        layout.addWidget(self.phi_k_text, 2, 5)

        layout.addWidget(self.C0_k_label, 3, 4)
        layout.addWidget(self.C0_k_text, 3, 5)

        layout.addWidget(self.pw_label, 3, 2)
        layout.addWidget(self.pw_text, 3, 3)

        layout.addWidget(self.z_label, 2, 0)
        layout.addWidget(self.z_text, 2, 1)

        layout.addWidget(self.ro_label, 3, 0)
        layout.addWidget(self.ro_text, 3, 1)

        layout.addWidget(self.m, 4, 1, 5, 5)
        # w.show()

        self.setLayout(layout)
        self.show()

        QtCore.QObject.connect(self.startButton, QtCore.SIGNAL("pressed()"), self.startAction)
        QtCore.QObject.connect(self.checkButton, QtCore.SIGNAL("pressed()"), self.checkAction)


    @QtCore.pyqtSignature("")
    def checkAction(self):
        logging.info("check started")
        self.startAction()

        # r_grid, theta_grid = np.meshgrid(r_vector, theta_vector)
        # gridShape = r_grid.shape

        # stress_at_well_list = self.window_calculate_section_at_well(777, self.angle_i_dec*np.pi/180, 0, self.pw,
        #                                                 self.sigma_H, self.sigma_H, self.sigma_v,
        #                                                 r_vector, theta_vector, gridShape,
        #                                                 Rw, self.z, 0.35)

        # stress_at_well_list = self.cur_section_cyl.stress_at_well
        # проверка критериями массива напряжений на стенке скважины
        # [True True True False False True False False]
        # False~Failure

        logging.info("general_stress_at_well_list:")
        coulomb_mohr1   = Criterion(RO_MIN, "FailureCoulombPrev")
        coulomb_mohr2   = Criterion(RO_MIN, "FailureCoulomb")
        coulomb_mohr3   = Criterion(RO_MIN, "CoulombMohr")
        drucker_prager1 = Criterion(RO_MIN, "DruckerPragerInnerCircle")
        drucker_prager2 = Criterion(RO_MIN, "DruckerPragerMiddleCircle")
        drucker_prager3 = Criterion(RO_MIN, "DruckerPragerOuterCircle")
        drucker_prager4 = Criterion(RO_MIN, "DruckerPragerPlainStrainState")
        drucker_prager5 = Criterion(RO_MIN, "DruckerPragerTriaxialCompression")

        criteria_list = []
        criteria_list.append(coulomb_mohr1)
        criteria_list.append(coulomb_mohr2)
        criteria_list.append(coulomb_mohr3)
        criteria_list.append(drucker_prager1)
        criteria_list.append(drucker_prager2)
        criteria_list.append(drucker_prager3)
        criteria_list.append(drucker_prager4)
        criteria_list.append(drucker_prager5)

        criteria_list_ans = [True for c in criteria_list]

        for stress_id, stress_cyl in enumerate(self.cur_section_cyl.stress_at_well):
            logging.info(str(stress_id)+": "+stress_cyl.getGeneralStressStr())
            logging.info("\ttrace comparison: "+str(stress_cyl.getTrace()-stress_cyl.get_gstress_trace()))
            logging.info("-------------------------------")
            logging.info(self.fcriteria.check_gstress(stress_cyl.getGeneralStress(), self.phi_k_dec*np.pi/180, self.C0_k))
            # logging.info("\t\t\tplot saved: "+str(stress_id)+".png")

            sigma3, sigma2, sigma1 = stress_cyl.getGeneralStress()
            plot_failure_criteria(stress_id, sigma3, sigma2, sigma1, self.phi_k_dec*np.pi/180, self.C0_k)
            logging.info("\t\t\tcoulomb_mohr1_ans="+    str(coulomb_mohr1.checkThisCriterion(stress_id,sigma3,sigma2,sigma1,self.phi_k_dec*np.pi/180,self.C0_k)))
            logging.info("\t\t\tcoulomb_mohr2_ans=" +   str(coulomb_mohr2.checkThisCriterion(stress_id, sigma3, sigma2, sigma1, self.phi_k_dec * np.pi / 180,self.C0_k)))
            logging.info("\t\t\tcoulomb_mohr3_ans=" +   str(coulomb_mohr3.checkThisCriterion(stress_id, sigma3, sigma2, sigma1, self.phi_k_dec * np.pi / 180,self.C0_k)))
            logging.info("\t\t\tdrucker_prager1_ans=" + str(drucker_prager1.checkThisCriterion(stress_id, sigma3, sigma2, sigma1, self.phi_k_dec * np.pi / 180,self.C0_k)))
            logging.info("\t\t\tdrucker_prager2_ans=" + str(drucker_prager2.checkThisCriterion(stress_id, sigma3, sigma2, sigma1, self.phi_k_dec * np.pi / 180,self.C0_k)))
            logging.info("\t\t\tdrucker_prager3_ans=" + str(drucker_prager3.checkThisCriterion(stress_id, sigma3, sigma2, sigma1, self.phi_k_dec * np.pi / 180,self.C0_k)))
            logging.info("\t\t\tdrucker_prager4_ans=" + str(drucker_prager4.checkThisCriterion(stress_id, sigma3, sigma2, sigma1, self.phi_k_dec * np.pi / 180,self.C0_k)))
            logging.info("\t\t\tdrucker_prager5_ans=" + str(drucker_prager5.checkThisCriterion(stress_id, sigma3, sigma2, sigma1, self.phi_k_dec * np.pi / 180,self.C0_k)))

            for id, criterion in enumerate(criteria_list):
                criteria_list_ans[id] = criteria_list_ans[id] and criterion.checkThisCriterion(stress_id,sigma3, sigma2, sigma1,self.phi_k_dec*np.pi/180, self.C0_k)

            logging.info("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            # logging.info("stress_cyl.getGeneralStressMtx()")
        # logging.info(stress_cyl.getGeneralStressMtx())
        logging.info("criterion_ans_at_section:")
        logging.info(criteria_list_ans)
        # fcriteria = FailureCriteria(RO_MAX, 777)
        # fcriteria_ans = fcriteria.check_section(self.cur_section_cyl.stress_at_well, self.ro, 777,
        #                                         self.phi_k_dec * np.pi / 180, self.C0_k)
        # logging.info(fcriteria_ans)


        # print(get_min_ro_in_section(777, self.angle_i_dec, 0, self.sigma_H, self.sigma_H, self.sigma_v,
        #                             r_vector, theta_vector, gridShape,
        #                             Rw, self.z, self.nu_k, self.phi_k, self.C0_k))


if __name__=='__main__':
    app = QtGui.QApplication(sys.argv)
    widget = MyWindow()
    widget.setGeometry(50,100, 1000,600)
    widget.show()
    sys.exit(app.exec_())
