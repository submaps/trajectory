import logging


class Layer:
    def __init__(self, id, z, h_k, ro_k, nu_k, phi_k, C0_k):
        self.id = int(id)
        self.z = float(z)
        self.h_k = float(h_k)
        self.ro_k = float(ro_k)
        self.nu_k = float(nu_k)
        self.phi_k = float(phi_k)
        self.C0_k = float(C0_k)

    def printLayer(self):
        logging.info(
            "layer " + str(self.id) + ": " + str(self.z) + " " + str(self.h_k) + " " + str(self.ro_k) + " " + str(
                self.nu_k) + " " + str(self.phi_k) + " " + str(self.C0_k))
