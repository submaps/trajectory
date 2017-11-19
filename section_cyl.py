import logging


class SectionCyl:
    def __init__(self, s_r, s_theta, s_z, t_r_theta, t_theta_z, t_r_z, i, a, pw, Rw, z, stress_at_well):
        self.s_r = s_r
        self.s_theta = s_theta
        self.s_z = s_z
        self.t_r_theta = t_r_theta
        self.t_theta_z = t_theta_z
        self.t_r_z = t_r_z
        # i, a, pw, Rw
        self.i = i
        self.a = a
        self.pw = pw
        self.Rw = Rw
        self.z = z
        self.stress_at_well = stress_at_well

    def printSectionCyl(self):
        logging.info("s_r       =" + str(self.s_r))
        logging.info("s_theta   =" + str(self.s_theta))
        logging.info("s_z       =" + str(self.s_z))
        logging.info("t_r_theta =" + str(self.t_r_theta))
        logging.info("t_theta_z =" + str(self.t_theta_z))
        logging.info("t_r_z     =" + str(self.t_r_z))
        logging.info("stress_at_well=" + str(self.stress_at_well))
