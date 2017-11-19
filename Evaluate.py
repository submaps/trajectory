# coding=utf-8
import numpy as np
import matplotlib.pyplot as pl
import random

def testEigvals():
     # test eigvals
     a = [[25000000, 0, 0],
          [0, 25000000, 0],
          [0, 0, 6000000]]
     print np.linalg.eigvals(a)
     print np.tan(np.pi / 4)

def generatePoints():
     # generate points
     l1 = 2000
     l2 = 1000
     l3 = 500
     l = 4 * l2
     r = l / (2 * np.pi)

     l3xstart = r
     l3ystart = l1 + r
     X1 = np.zeros(20)
     Y1 = np.linspace(0, -l1, 20)

     X3 = np.linspace(int(l3xstart), int(l3xstart + l3), 20)
     Y3 = np.linspace(-l3ystart, -l3ystart, 20)

     indent = 2000
     x0 = r
     y0 = r

     print "r=", r
     X2 = np.linspace(1, int(r), int(r / 25))
     # print X
     Y2 = map(lambda x: y0 - (r ** 2 - (x - x0) ** 2) ** 0.5 - l1 - r, X2)

     fig = pl.figure()
     ax = fig.add_subplot(111, aspect='equal')
     ax.plot(X1, Y1, color='blue', marker="^")
     ax.plot(X2, Y2)
     ax.plot(X3.flat, Y3.flat, color="red", marker="o")
     pl.show()

     for (x, y) in zip(X1, Y1):
          print(str(x) + ";" + str(y) + ";")

     for (x, y) in zip(X2, Y2):
          print(str(x) + ";" + str(y) + ";")

     for (x, y) in zip(X3.flat, Y3.flat):
          print(str(x) + ";" + str(y) + ";")

def testBoolAnd():
     criteria_list_id = range(8)
     # criteria_list_ans = [True for c in criteria_list]
     # criteria_list_ans[2]=False
     #
     # criteria_list_ans=[True]*len(criteria_list_id)
     # # создаем копии
     # for criterion_id in criteria_list_id:
     #      cur_ans = True
     #      if criterion_id%2 == 0:
     #           cur_ans = False
     #
     #           # print False
     #      # criterion_ans = bool(criterion_ans) and bool(cur_ans)
     #      criteria_list_ans[criterion_id] = criteria_list_ans[criterion_id] and cur_ans
     #
     # print criteria_list_ans

     for id in criteria_list_id:
          if id%2==0: id=2

     print criteria_list_id

     # criteria_list = range(8)
     # criteria_list_ans = [1 for c in criteria_list]
     #
     # for criterion_id, criterion_ans in enumerate(criteria_list_ans):
     #      cur_ans=1
     #      if criterion_id%2 == 0: cur_ans=0
     #      criterion_ans *= cur_ans
     #      print str(criterion_id) + ": " + str(criterion_ans)
     #
     # print criteria_list_ans

     # 0: 0
     # 1: 1
     # 2: 0
     # 3: 1
     # 4: 0
     # 5: 1
     # 6: 0
     # 7: 1
     # [1, 1, 1, 1, 1, 1, 1, 1]

def testRoGraph():
     # -*- coding: UTF-8 -*-
     import numpy
     import pylab
     import matplotlib.ticker
     xvals = numpy.arange(-10.0, 10.1, 0.1)
     yvals = numpy.sinc(xvals)
     figure = pylab.figure()
     axes = figure.add_subplot(1, 1, 1)
     # Создаем форматер
     formatter = matplotlib.ticker.NullFormatter()
     # Установка форматера для оси X
     axes.xaxis.set_major_formatter(formatter)
     pylab.plot(xvals, yvals)
     axes.grid()
     pylab.show()


     # ax = fig.add_subplot(111)



# testBoolAnd()
testRoGraph()