#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '08/27/2021'
__version__ = '1.0'




PN = 0.012   # natural abundance of 13C


import numpy as np
from functools import lru_cache




class NMRrun:
	
	def __init__(self, P1):
		'''
		Parameters
		P1: float
			the determined overall degree of 13C labeling
		'''
		
		self.P1 = P1
		self.Pn = PN
		
		self.Pf = (self.P1 - self.Pn)/(1 - self.Pn)
	
	
	@property
	@lru_cache()
	def K_termC(self):
		
		K_1_s = 1 - self.P1
		K_1_d = self.P1
		
		K_2_s = self.Pn*(1 - self.Pn)*(1 - self.Pf)/self.P1
		K_2_d = (self.Pn**2*(1 - self.Pf) + self.Pf)/self.P1
		
		return np.array([[K_1_s, K_2_s], 
						 [K_1_d, K_2_d]])
	
	
	@property
	@lru_cache()	
	def K_centC_ne(self):
		
		K_1_s = (1 - self.P1)**2
		K_1_da = self.P1*(1 - self.P1)
		K_1_db = K_1_da
		K_1_dd = self.P1**2
		
		K_2a_s = self.Pn*(1 - self.Pn)*(1 - self.P1)*(1 - self.Pf)/self.P1
		K_2a_da = (self.Pn**2*(1 - self.P1)*(1 - self.Pf) + self.Pf*(1 - self.P1))/self.P1
		K_2a_db = self.Pn*(1 - self.Pn)*(1 - self.Pf)
		K_2a_dd = self.Pn**2*(1 - self.Pf) + self.Pf
		
		K_2b_s = K_2a_s
		K_2b_da = K_2a_db
		K_2b_db = K_2a_da
		K_2b_dd = K_2a_dd
		
		K_3_s = self.Pn*(1 - self.Pn)**2*(1 - self.Pf)/self.P1
		K_3_da = self.Pn**2*(1 - self.Pn)*(1 - self.Pf)/self.P1
		K_3_db = K_3_da
		K_3_dd = (self.Pn**3*(1 - self.Pf) + self.Pf)/self.P1
		
		return np.array([[K_1_s, K_2a_s, K_2b_s, K_3_s],
						 [K_1_da, K_2a_da, K_2b_da, K_3_da],
						 [K_1_db, K_2a_db, K_2b_db, K_3_db],
						 [K_1_dd, K_2a_dd, K_2b_dd, K_3_dd]])
		
	
	@property
	@lru_cache()
	def K_centC_eq(self):
		
		K_1_s = (1 - self.P1)**2
		K_1_d = self.P1*(1 - self.P1)*2
		K_1_t = self.P1**2
		
		K_2_s = self.Pn*(1 - self.Pn)*(1 - self.P1)*(1 - self.Pf)/self.P1*2
		K_2_d = (self.Pn**2*(1 - self.P1)*(1 - self.Pf) + self.Pf*(1 - self.P1))/self.P1 + self.Pn*(1 - self.Pn)*(1 - self.Pf)
		K_2_t = (self.Pn**2*(1 - self.Pf) + self.Pf)*2
		
		K_3_s = self.Pn*(1 - self.Pn)**2*(1 - self.Pf)/self.P1
		K_3_d = self.Pn**2*(1 - self.Pn)*(1 - self.Pf)/self.P1*2
		K_3_t = (self.Pn**3*(1 - self.Pf) + self.Pf)/self.P1
		
		return np.array([[K_1_s, K_2_s, K_3_s],
						 [K_1_d, K_2_d, K_3_d],
						 [K_1_t, K_2_t, K_3_t]])
		




