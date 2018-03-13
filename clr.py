import glob
import numpy as np
import rafias_lib as rl
import matplotlib.pyplot as plt

files_clr1 = np.sort(glob.glob('/data1/tso_analysis/all_tso_cv3/raw_separated_MMM/NRCN821CLRSUB1*_481_SE_2016-*/*.slp.fits'))

files_clr2 = np.sort(glob.glob('/data1/tso_analysis/all_tso_cv3/raw_separated_MMM/NRCN821CLRSUB1*_489_SE_2016-*/*.slp.fits'))

centers_clr1 = rl.gen_center_g2d(164,161,5,3500,2,2,0, files_clr1)

centers_clr2 = rl.gen_center_g2d(160,155,5,4500,2,2,0, files_clr2)

rt_clr = rl.radius_testing(centers_clr1[1], centers_clr1[1], files_clr1, centers_clr2[1], centers_clr2[1], files_clr2, 1, 5, 1, 5, 10, 1, 10, 15, 1)

ascii.write(rt_clr, 'rt_clr.csv', overwrite = True) 
