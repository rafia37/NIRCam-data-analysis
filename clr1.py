import glob
import numpy as np
import rafias_lib as rl
import matplotlib.pyplot as plt

r1, rin1, rout1 = raw_input('r,rin,rout: ').split(',')
r, rin, rout = int(r1), int(rin1), int(rout1)
print 'r = %i, rin = %i, rout = %i' % (r, rin, rout)

#"Analyzing In-focus data" Notebook
 
a1_files = np.sort(glob.glob('/data1/tso_analysis/all_tso_cv3/raw_separated_MMM'
                              '/NRCN821CLRSUB1*_481_SE_2016-*/*.slp.fits'))

b4_files = np.sort(glob.glob('/data1/tso_analysis/all_tso_cv3/raw_separated_MMM'
                              '/NRCN821CLRSUB1*_489_SE_2016-*/*.slp.fits'))

a1_centers = rl.gen_center_g2d(164,161,5,3500,2,2,0, a1_files)

b4_centers = rl.gen_center_g2d(160,155,5,4500,2,2,0, b4_files)

print 'generated 1st set centers'

a1_data = rl.time_series(a1_centers[1], a1_centers[2], a1_files, r, rin, rout)

b4_data = rl.time_series(b4_centers[1], b4_centers[2], b4_files, r, rin, rout)

print 'generated 1st set time series'

av_data = (a1_data['res_flux']+b4_data['res_flux'])/2 

st_dev = np.std(av_data[3:]/np.median(av_data[3:]))

print 'original stdev:', st_dev


#"Radius Testing notebook"

files_clr1 = np.sort(glob.glob('/data1/tso_analysis/all_tso_cv3/'
                               'raw_separated_MMM/NRCN821CLRSUB1*_481_SE_2016-*/*.slp.fits'))

files_clr2 = np.sort(glob.glob('/data1/tso_analysis/all_tso_cv3/'
                               'raw_separated_MMM/NRCN821CLRSUB1*_489_SE_2016-*/*.slp.fits'))

centers_clr1 = rl.gen_center_g2d(164,161,5,3500,2,2,0, files_clr1)

centers_clr2 = rl.gen_center_g2d(160,155,5,4500,2,2,0, files_clr2)

print 'generated 2nd set centers'

data1 = rl.time_series(centers_clr1[1], centers_clr1[2], files_clr1, r, rin, rout)

data2 = rl.time_series(centers_clr2[1], centers_clr2[2], files_clr2, r, rin, rout)

print 'generated 2nd set time series'

detrended1 = rl.linear_bestfit(data1['time'], data1['res_flux'], 0.00002, 1)

detrended2 = rl.linear_bestfit(data2['time'], data2['res_flux'], 0.00002, 1)

av = (detrended1 + detrended2)/2.

stdev = np.std(av/np.median(av))

print 'rad_test stdev:', stdev

plt.plot(a1_data['time'], av_data/np.median(av_data), label = 'original: %f ppm' % (st_dev*1000000))
plt.plot(data1['time'], av/np.median(av), label = 'rad test: %f ppm' % (stdev*1000000))
plt.legend(loc = 'best')
plt.xlabel('time[sec]')
plt.ylabel('flux [DN/s]')
plt.title('%i, %i, %i' % (r, rin, rout))
plt.show()


