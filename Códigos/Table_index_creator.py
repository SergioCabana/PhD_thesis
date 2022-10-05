# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 12:46:26 2022

@author: caban

objetivo: array de numpy con indices de tablas de AIRES 19.04.08
las tablas con 9999 no existen
"""

import numpy as np

tables = np.ones((23,16), dtype = int)

col0 = np.array([1001, 1005, 1006, 1007, 1008, 1011, 1012, 1013, 1014, 1021, \
                 1022, 1023, 1041, 1091, 1092, 1205, 1207, 1211, 1213, 1291, \
                 1292, 1293, 1296], dtype = int).reshape((23,))
    
tables[:,0] = col0

tables[:15, 1] = tables[:15,0]+300
tables[15:, 1] = tables[15:,0]+200

tables[:,2] = tables[:,0]+500

tables[:,3] = tables[:,0]+1000

tables[:15, 4] = tables[:15,0]+1300
tables[15:, 4] = tables[15:,0]+1200

tables[:,5] = tables[:,0]+1500

tables[:15, 6] = tables[:15,0]+1800
tables[15:, 6] = tables[15:,0]+1700

tables[:,7] = tables[:,0]+2000

tables[:,8] = tables[:,0]+4000

tables[:,9] = tables[:,0]+5000

tables[:15, 10] = tables[:15,0]+5300
tables[15:, 10] = tables[15:,0]+5200

tables[:,11] = tables[:,0]+5500

tables[:,12] = tables[:,0]+6000

tables[:15, 13] = tables[:15,0]+6300
tables[15:, 13] = tables[15:,0]+6200

tables[:,14] = tables[:,0]+6500

tables[:15, 15] = tables[:15,0]+6800
tables[15:, 15] = tables[15:,0]+6700

tables[22,:9] = tables[22,12:] = 9999
tables[1:14,7] = tables[17:19,7] = 9999
tables[5:13,12:] = 9999
tables[17:19,12:] = 9999

np.savetxt('AIRES_table_index.txt', tables, fmt = '%4d')