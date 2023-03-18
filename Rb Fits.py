import matplotlib.pyplot as plt
import numpy as np
import csv
from glob import glob
import peakutils
from scipy.optimize import curve_fit
from scipy.special import wofz,erf
from lmfit.models import Model
import time

start=time.time()




plt.style.use('default')

filenames=glob('spectrum'+'*.csv')
# filenames=glob('*.csv')
# for i in range(10):
	# print(filename[50+i])

# filenames=filenames[83:]
# filenames=filenames[0:3]

def voigt(xdata,amplitude,center,sigma,gamma,offset):
	# sigma=np.abs(sigma)
	# gamma=np.abs(gamma)
	z=[(val-center+1j*gamma)/(sigma*np.sqrt(2)) for val in xdata]
	w=wofz(z)
	f=(-amplitude*np.real(w))/(sigma*np.sqrt(2*np.pi))+offset
	return f

def linear_function(xdata,m,b):
	return [m*num+b for num in xdata]

def gaussian(xdata,amplitude,center,variance,offset):
	f=[-amplitude*np.exp(-0.5*((val-center)/variance)**2)+offset for val in xdata]
	return f

def remove_slope(xdata,ydata):
	slope=(ydata[-1]-ydata[0])/(xdata[-1]-xdata[0])
	slope_vals=[slope*elem for elem in xdata]
	yvals=np.subtract(ydata,slope_vals)
	return yvals

def normalize(array):
	new_array=[]
	maximum=max(array)
	minimum=min(array)
	for i in array:
 		new_array.append((i)/(maximum))
	return new_array

def convert_to_frequencies(xdata,ydata,frequency_calib_list):
	ydata=np.asarray(ydata)	
	indexes = peakutils.indexes(ydata, thres=0.2, min_dist=40)
	print(indexes)
	frequencies=[]
	times=[xdata[index] for index in indexes]
	if len(indexes)==len(frequency_calib_list):
		fit=np.poly1d(np.polyfit(times,frequency_calib_list,1))
		frequencies=fit(xdata)
	ydata=[-val for val in ydata]
	frequencies=[-val for val in frequencies]	
	return frequencies


def sum_of_voigts(xdata,a1,c1,a2,c2,s,v,offset):
	sum=[a+b for a,b in zip(voigt(xdata,a1,c1,s,v,offset),voigt(xdata,a2,c2,s,v,0))]
	# sum=gaussian(xdata,a1,c1,v,offset)+gaussian(xdata,a2,c2,v,0)
	# print(sum)
	return sum

def voigt_linewidth(sigma,gamma):
	#sigma is the lorentzian linewidth
	#gamma is the gaussian variance
	fg=2*gamma*np.sqrt(2*np.log(2))
	fl=2*sigma
	fv=fl/2+np.sqrt(fl**2/4+fg**2)
	return fv

# def calibrate():
# 	plt.plot(x,ysig)
# 	plt.plot(x,yref)
# 	#plt.plot(x,fpcavity)
# 	plt.xlim(min(x),max(x))
# 	plt.ylim(0,1)
# 	print("please click on the Rb 85 and Rb87 peaks")
# 	points=plt.ginput(2)
# 	plt.close()
# 	return points

def calibrate():
	axs[0].plot(x,ysig)
	axs[0].plot(x,yref)
	#plt.plot(x,fpcavity)
	# axs[0].xlim(min(x),max(x))
	# axs[0].ylim(0,1)
	print("please click on the Rb 85 and Rb87 peaks")
	points=plt.ginput(2)
	plt.close()
	return points

def smoother(array,N):
	new_array=np.convolve(array,np.ones(N)/N)
	new_array=new_array[N:-N]
	return new_array

def sloped_voigt(xdata,amplitude,center,sigma,gamma,offset,m,b):
	new_array = [sum(i) for i in zip(voigt(xdata,amplitude,center,sigma,gamma,offset),linear_function(xdata,m,b))] 
	return new_array

def sum_of_Rb_voigts(xdata,a1,c1,a2,c2,s,v,m):
	#Rb85 is 72.17% abundant
	#Rb87 is 27.83% abundant
	#Relative amplitudes should be the ratio of these factors

	# rel_ratio=0.2783/0.7217
	rel_ratio=1


	#Rb85 and Rb87 hav doppler broadened resonances:
	#From Siddons 2008
	#Isotopic degeneracy yields 1/12 for 85 and 1/8 for 87
	#F=3 to F'=4 --> 0 MHz
	#F=3 to F'=3 --> -121 MHz
	#F=3 to F'=2 --> -184 MHz
	Rb85peak3to4 = voigt(xdata,a1*1*(1/12),c1-0,np.abs(s),np.abs(v),0)
	Rb85peak3to3 = voigt(xdata,a1*35/81*(1/12),c1-0.121,np.abs(s),np.abs(v),0)
	Rb85peak3to2 = voigt(xdata,a1*10/81*(1/12),c1-0.184,np.abs(s),np.abs(v),0)
	#F=2 to F'=3 --> 0 MHz
	#F=2 to F'=2 --> -267 MHz
	#F=2 to F'=1 --> -424 MHz
	Rb87peak2to3 = voigt(xdata,a2*(7/9)*rel_ratio*(1/8),-1.1-0+c2,np.abs(s),np.abs(v),0)
	Rb87peak2to2 = voigt(xdata,a2*(5/18)*rel_ratio*(1/8),-1.1-0.267+c2,np.abs(s),np.abs(v),0)
	Rb87peak2to1 = voigt(xdata,a2*(1/18)*rel_ratio*(1/8),-1.1-0.424+c2,np.abs(s),np.abs(v),0)

	log_laser_slope=[np.log(num) for num in linear_function(xdata,m,1)]

	new_array=np.array([Rb85peak3to4,Rb85peak3to3,Rb85peak3to2,Rb87peak2to3,Rb87peak2to2,Rb87peak2to1,log_laser_slope])
	# print(np.shape(new_array.sum(axis=0)))
	return new_array.sum(axis=0)

def sum_of_Rb_Gaussians(xdata,a1,a2,c1,v,m):

	#Rb85 and Rb87 have doppler broadened resonances:
	#From Siddons 2008

	Rb85peakFg3 = gaussian(xdata,a1*3,c1-0,np.abs(v),0)

	# Isotopic separation from line center is 1.1 GHz

	Rb87peakFg2 = gaussian(xdata,a2*2,-1.1-0+c1,np.abs(v),0)

	log_laser_slope=[np.log(num) for num in linear_function(xdata,m,1)]

	new_array=np.array([Rb85peakFg3,Rb87peakFg2,log_laser_slope])
	# print(np.shape(new_array.sum(axis=0)))
	return new_array.sum(axis=0)

linewidths=[]
shifts=[]
voigt_measured_linewidths=[]

fits=[]



fig, axs=plt.subplots(2)
# fig.suptitle('Spectra and their linewidths \n 90 C and 75 microW')

initial_calibration=1

for file in filenames:
	# with open(str(text)+'_wld_.csv','r') as f:
	x=[]
	ysig=[]
	peaks=[]
	yref=[]
	print('\n')
	print('Analyzing file '+str(file))
	print('\n')
	with open(file,'r') as f:
		reader=csv.reader(f)
		i=0
		for row in reader:
			if i>5:
				row=row[0].split(' ')
				# print(row)
				# print(row)
				# if i>5:
				x.append(float(row[0]))
				ysig.append(float(row[1]))
				yref.append(float(row[3]))
			i=i+1
		# ysig=remove_slope(x,ysig)
		# yref=remove_slope(x,yref)

		freq_points=[[0.27,2.4],[2.4,0.27]]
		# freq_points=[[0.17,1.5],[2.23,3.29]]
		freq_points=[[2.23,2.23],[0.16,0.16]]
		if initial_calibration==0:
			freq_points=calibrate()
			initial_calibration=1
			print(freq_points)
		x=[(1.1/(freq_points[1][0]-freq_points[0][0]))*(num-freq_points[1][0]) for num in x]

		#Normalize for consistency
		ysig=normalize(ysig)
		yref=normalize(yref)

		#Take ln to account for Beer's Law
		yref=[np.log(num/np.max(yref)) for num in yref]
		ysig=[np.log(num/np.max(ysig)) for num in ysig]		

		sObject=slice(700,1800)	

		# Smooth data if necessary
		# moving_avg_size=1
		# yref=np.convolve(yref,np.ones(moving_avg_size)/moving_avg_size)
		# ysig=np.convolve(ysig,np.ones(moving_avg_size)/moving_avg_size)
		# x=np.convolve(x,np.ones(moving_avg_size)/moving_avg_size)

		# Slice data if necessary
		# yref=yref[sObject]
		# ysig=ysig[sObject]
		# x=x[sObject]

		# ysig=yref


		gmodel=Model(sum_of_Rb_Gaussians)

		params=gmodel.make_params()

		params.add('a1',min=0.1,max=20)
		params.add('c1',min=-0.5,max=0.5)			
		params.add('a2',min=0.1,max=10)
		# params.add('c2',min=-0.2,max=0.2)
		# params.add('s',min=0.01,max=1)
		params.add('v',min=0.2,max=3)
		params.add('m',min=-0.2,max=0.2)
		# params.add('offset',min=0,max=1)

		# result=gmodel.fit(ysig,params,xdata=x,a1=5,c1=0,a2=3,c2=0,s=1,v=0.2,m=0,fit_kws={'maxfev': 10})
		# result=gmodel.fit(ysig,params,xdata=x,a1=20,c1=0,a2=7,c2=0,s=0.01,v=0.3,m=0,fit_kws={'xtol': 1.e-7,'ftol':1.e-7,'maxfev':1000})
		# result=gmodel.fit(ysig,params,xdata=x,a1=15,c1=0,a2=7,c2=0,s=0.01,v=0.3,m=0,fit_kws={'maxfev':2000})
		result=gmodel.fit(ysig,params,xdata=x,a1=2,a2=1,c1=0,v=0.3,m=0,fit_kws={'maxfev':2000})
		# result=gmodel.fit(ysig,params,xdata=x,a1=5,c1=0,a2=3,c2=0,s=1,v=0.2,m=0)

		# axs[0].plot(x,result.init_fit)
		# axs[0].plot(x,result.best_fit,'r:')
		# axs[0].plot(x,result.best_fit,alpha=0.5,color='#000000')

		fits.append([x,result.best_fit])

		# np.savetxt('best_fit_of_'+file,np.transpose([x,result.best_fit]),delimiter=',')

		# plt.plot(x,result.best_fit,'r:')

		coeff=result.params

		if coeff['a1']==0.1 and coeff['a2']==0.1:
			print('No Rb in this one')

		print('Linewidth is '+str(2*np.sqrt(2*np.log(2))*np.round(coeff['v'].value,2)))
		# shifts.append()
		voigt_measured_linewidths.append(2*np.sqrt(2*np.log(2))*coeff['v'])
		result.params.pretty_print()

	# axs[0].plot(x,yref)
	axs[0].plot(x,ysig,color='#000000',alpha=0.5)
	# axs[0].plot(x,ysig,color='#000000')
	# plt.plot(x,ysig)

	# axs[0].plot(x,yref)

num_cells_with_rb=0
cells_with_rb=[]

for line in linewidths:
	if line>0.001:
		num_cells_with_rb+=1
		cells_with_rb+=line



def center_spectra(array,min_index):
	new_array=[]
	for element in array:
		new_array.append(element-array[int(min_index)])
	return new_array

# axs[0].set_title('BWLD1 Optical Cavity Spectroscopy \n 75 microW at 90 C')
# plt.ylabel('Normalized transmission')
# plt.xlabel('Laser frequency (GHz)')

axs[0].set_title('BWLD1 Optical Cavity Spectroscopy \n 1 mm beam with 75 microW at 90 C \n Fit Using Two Gaussians')
axs[0].set_title('Optical Spectroscopy of Vapor Cells')
axs[0].set_xlabel('Laser Frequency Detuning (GHz)')
axs[0].set_ylabel('Logarithm of \n Normalized Transmission')


# plt.title('BWLD1 Optical Cavity Spectroscopy \n 75 microW at 90 C \n Gaussian Fits Hi Res')
# plt.xlabel('Laser Frequency Detuning (GHz)')
# plt.ylabel('Normalized Transmission')

bin_num=20
Torr_tot=100

# axs[1].set_title('Histogram of Buffer-Gas-Enhanced Linewidths with '+str(bin_num)+' bins')
# axs[1].set_xlabel('Linewidth Minus Doppler Broadening (GHz)')
# axs[1].set_ylabel('Occurences')
# axs[1].set_xticks(np.linspace(0,1,9))



# axs[1].set_title('Histogram of Residual Pressures in Cells with '+str(bin_num)+' bins')
axs[1].set_title('Histogram of Residual Pressures in Cells ('+str(bin_num)+' bins)')
axs[1].set_xlabel('Pressure (Torr)')
axs[1].set_ylabel('Occurences')
axs[1].set_xticks(np.linspace(0,Torr_tot,bin_num+1))

fig.subplots_adjust(hspace=1)

# from rotondaro 1997, the shift is -5.79 MHz/Torr 
# and the broadening is 18.3 MHz/Torr for the Rb D2 transition

linewidths=voigt_measured_linewidths

broadening=[val-min(linewidths) for val in linewidths]

pressures=[num*1000/18.3 for num in broadening]

np.savetxt('linewidths.csv',linewidths)

np.savetxt('pressures.csv',pressures)

# np.savetxt('fits.csv',fits,delimiter=',')


print('this script took '+str(time.time()-start)+' s to run '+str(len(filenames))+' fits')


# print(len(filenames)-len(linewidths))
# axs[1].hist(broadening,bins=np.linspace(0,0.5,bin_num),color='#EE6666')
# axs[1].hist(pressures,bins=np.linspace(0,25,bin_num),color='#66ee66')
axs[1].hist(pressures,bins=np.linspace(0,Torr_tot,bin_num),color='#000000',alpha=0.7,rwidth=0.9)
plt.savefig('WLD13_residual_pressures.png')
plt.show()
