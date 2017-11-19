import logging


class Point:
    def __init__(self, id, x, y, z):
        self.id = float(id)
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def printPoint(self):
        logging.info("point " + str(self.id) + ": " + str(self.x) + " " + str(self.y) + " " + str(self.z))

