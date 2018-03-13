from astropy.io import ascii
import numpy as np
import glob, pdb
import rafias_lib as rl


files_sub1 = np.sort(glob.glob('/data1/tso_analysis/all_tso_cv3/raw_separated_MMM/NRCN821WLP8SUB-60*_1_481_SE_*/*.slp.fits'))
files_sub2 = np.sort(glob.glob('/data1/tso_analysis/all_tso_cv3/raw_separated_MMM/NRCN821WLP8SUB-60*_1_489_SE_*/*.slp.fits'))[:306]

print "files loaded"
centers_sub1 = rl.gen_center_g2d(166,160,5,400,2,2,0,files_sub1)
print "Done with centers_sub1"
centers_sub2 = rl.gen_center_g2d(162,156,5,500,2,2,0,files_sub2)

print "Entering radius testing loop sub320"

rt_sub = rl.radius_testing(centers_sub1[1], centers_sub1[2], files_sub1, centers_sub2[1], centers_sub2[2], files_sub2, 60, 80, 2, 70, 90, 2, 80, 100, 2)

print "Done with sub320"

rt_sub.write('rt_sub(fg).csv')





files_sub6401 = np.sort(glob.glob('/data1/tso_analysis/all_tso_cv3/raw_separated_MMM/NRCN821WLP8SUB640-60*_1_481_SE_*/*.slp.fits'))
files_sub6402 = np.sort(glob.glob('/data1/tso_analysis/all_tso_cv3/raw_separated_MMM/NRCN821WLP8SUB640-60*_1_489_SE_*/*.slp.fits'))

print "files loaded sub640"
centers_sub6401 = rl.gen_center_g2d(326,320,5,400,2,2,0,files_sub6401)
print "Done with centers_sub6401"
centers_sub6402 = rl.gen_center_g2d(319,316,5,500,2,2,0,files_sub6402)

print "Entering radius testing loop sub640"

rt_sub640 = rl.radius_testing(centers_sub6401[1], centers_sub6401[2], files_sub6401, centers_sub6402[1], centers_sub6402[2], files_sub6402, 50, 70, 2, 60, 80, 2, 70, 90, 2) 

print "Done with sub640"

rt_sub640.write('rt_sub640(fg).csv')






files_full1 = np.sort(glob.glob('/data1/tso_analysis/all_tso_cv3/raw_separated_MMM/NRCN821WLP8FULL1-*_1_481_SE_*/*.red.fits'))
files_full2 = np.sort(glob.glob('/data1/tso_analysis/all_tso_cv3/raw_separated_MMM/NRCN821WLP8FULL1-*_1_489_SE_*/*.red.fits'))

print "files loaded full"
centers_full1 = rl.gen_center_g2d(1405,1036,5,400,4,4,0, files_full1, red = True)
print "Done with centers_full"
centers_full2 = rl.gen_center_g2d(828,821,5,600,4,4,0, files_full2, red = True)

print "Entering radius testing loop full"

rt_full = rl.radius_testing(centers_full1[1], centers_full1[2], files_full1, centers_full2[1], centers_full2[2], files_full2, 30, 90, 10, 50, 110, 10, 70, 130, 10, Red = True)  

print "Done with full"

rt_full.write('rt_full.csv')







files_clr1 = np.sort(glob.glob('/data1/tso_analysis/all_tso_cv3/raw_separated_MMM/NRCN821CLRSUB1*_481_SE_2016-*/*.slp.fits'))
files_clr2 = np.sort(glob.glob('/data1/tso_analysis/all_tso_cv3/raw_separated_MMM/NRCN821CLRSUB1*_489_SE_2016-*/*.slp.fits'))

print "files loaded clr"
centers_clr1 = rl.gen_center_g2d(164,161,5,3500,2,2,0, files_clr1)
print "Done with centers_clr"
centers_clr2 = rl.gen_center_g2d(160,155,5,4500,2,2,0, files_clr2)

print "Entering radius testing loop clr"

rt_clr = rl.radius_testing(centers_clr1[1], centers_clr1[1], files_clr1, centers_clr2[1], centers_clr2[1], files_clr2, 2, 6, 0.2, 4, 8, 0.2, 6, 10, 0.2)

print "Done with clr"   

rt_clr.write('rt_clr(fg).csv')
