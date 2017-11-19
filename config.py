# -*- coding: utf-8 -*-
#project config
import numpy as np
import logging
#уровень вывода лога

# log_level = logging.INFO
log_level = logging.DEBUG

# входные данные
input_trajectory = "experiment/trajectory.csv"
input_layers = "experiment/layers.csv"

# папка для выходных данных
output_dir = "experiment"

# density of drilling fluid kg/m3
# ro_analysis()
RO_MIN = 800
RO_MAX = 1600

# density of drilling fluid kg/m3
# calculate_each_sections()
ro_init = 1000

g = 9.81
Rw = 0.01


# sigma_v = 0
# sigma_H = 0

theta_vector = np.linspace(0, 2 * np.pi, 100)
r_vector = np.linspace(Rw, Rw * 2, 25)

target = "ro_analysis"
# target="calculate_trajectory"