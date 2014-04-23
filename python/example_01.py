#!/usr/bin/env python

import Image
import scipy.signal as ss
import numpy as np
import random
import matplotlib.pyplot as plt

def synthesize(imin,sizeout,tilesize,overlap):

	imout = np.zeros(sizeout,dtype=np.float64)
	sizein = imin.shape


	temp = np.ones((overlap,tilesize))
	errtop = ss.correlate2d(imin**2,temp)
	temp = np.ones((tilesize,overlap))
	errside = ss.correlate2d(imin**2,temp)
	
	temp = np.ones((tilesize-overlap,overlap))
	errsmall = ss.correlate2d(imin**2,temp)
	
	scale = 1.21
	
	for i in range(0,sizeout[0]-tilesize+1,tilesize-overlap):
		for j in range(0,sizeout[1]-tilesize+1,tilesize-overlap):
			
			if i>0 and j>0:
				shared = imout[i:i+overlap,j:j+tilesize]
				err = errtop -2*ss.correlate2d(imin,shared)+sum(sum(shared**2))
				errsize=err.shape
				err = err[overlap-1:errsize[0]-tilesize,tilesize-1:errsize[1]-tilesize]

				shared = imout[i+overlap:i+tilesize,j:j+overlap]
				err2 = errsmall-2*ss.correlate2d(imin,shared)+sum(sum(shared**2))
				err2size = err2.shape
				err = err+ err2[tilesize-1:err2size[0]-tilesize+overlap,overlap-1:err2size[1]-tilesize]

				mindata = err.min()
				threshold = scale*mindata
				
				best = np.nonzero(err<=threshold)
				c = random.randint(0,len(best[0])-1)
				pos = [best[0][c],best[1][c]]

			elif i>0:
				shared = imout[i:i+overlap,j:j+tilesize]
				err = errtop - 2*ss.correlate2d(imin,shared)+sum(sum(shared**2))
				errsize = err.shape
				err = err[overlap-1:errsize[0]-tilesize+1,tilesize-1:errsize[1]-tilesize+1]
				best = []
				mindata = err.min()
				threshold = scale*mindata
				
				best = np.nonzero(err<=threshold)
				c = random.randint(0,len(best[0])-1)
				pos = [best[0][c],best[1][c]]
			elif j>0:
				shared = imout[i:i+tilesize,j:j+overlap]
				err = errside - 2*ss.correlate2d(imin,shared)+sum(sum(shared**2))
				errsize = err.shape
				err = err[tilesize-1:errsize[0]-tilesize+1,overlap-1:errsize[1]-tilesize+1]
				mindata = err.min()
				threshold = scale*mindata
				
				best = np.nonzero(err<=threshold)
				c = random.randint(0,len(best[0])-1)
				pos = [best[0][c],best[1][c]]
			else:
				pos = [random.randint(0,sizein[0]-tilesize),random.randint(0,sizein[1]-tilesize)]
				
			
			for indexi in range(tilesize):
				for indexj in range(tilesize):
					imout[i+indexi,j+indexj] = imin[pos[0]+indexi,pos[1]+indexj]


	return imout


def main():
	img =Image.open('D.gif')
	array = np.array(img,dtype=float)
	sizeout = [array.shape[0]*4,array.shape[1]*4]
	
	imout=synthesize(array,sizeout,12,2)
	imout = imout.astype(np.uint8)
	

	img.show()
	result=Image.frombuffer('L',imout.shape,imout,'raw','L',0,1)
	result.show()

if __name__ == '__main__':
	main()


