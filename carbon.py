#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '08/25/2021'
__version__ = '1.0'




import numpy as np
import pandas as pd
from scipy.linalg import solve




class BaseC:
	
	def __init__(self, name, I):
		'''
		Parameters
		name: str
			name of carbon atom, e.g. "Ala_alpha"
		I: array or list of float
			observed relative multiplet intensities
		'''
		
		self.I = np.array(I)
		self.I[self.I < 0] = 0
		self.I = self.I/self.I.sum()
		
		
	def solve_f(self, trans_matrix, standardize = True):
		'''
		Parameters
		trans_matrix: array
			transforamtion matrix K
		'''
		
		self.F = solve(trans_matrix, self.I)
		
		if standardize:
			self.F[self.F < 0] = 0
			self.F = self.F/self.F.sum()
		
		return self.F


class TermC(BaseC):
	'''
	Terminal carbon atom with a singlet "s" and a doublet "d"
	'''
	
	def __init__(self, name, I):
		'''
		Parameters
		name: str
			name of carbon atom, e.g. "Ala_alpha"
		I: array or list of float
			observed relative multiplet intensities of "s" and "d"
		'''
		
		super().__init__(self, I)
		
		if self.I.size != 2:
			raise ValueError('The length of multiplet intensities of terminal carbon atom should be 2')
			
	
	@property
	def f1(self):
		
		if hasattr(self, 'F'):
			return round(self.F[0], 4)
			
		else:
			raise AttributeError('run solve_f first')
	
	
	@property
	def f2(self):
		
		if hasattr(self, 'F'):
			return round(self.F[1], 4)
			
		else:
			raise AttributeError('run solve_f first')
			

class CentCne(BaseC):
	'''
	Central carbon atom with a singlet "s", a doublet with small coupling constant "da",
	a doublet with big coupling constant "db", and a doublet of doublets "dd"
	'''
	
	def __init__(self, name, I):
		'''
		Parameters
		name: str
			name of carbon atom, e.g. "Ala_alpha"
		I: array or list of float
			observed relative multiplet intensities of "s", "da", "db" and "dd"
		'''
		
		super().__init__(self, I)
		
		if self.I.size != 4:
			raise ValueError('The length of multiplet intensities of central carbon atom with unequal coupling constants to the adjacent carbons should be 4')
		
	
	@property
	def f1(self):
		
		if hasattr(self, 'F'):
			return round(self.F[0], 4)
			
		else:
			raise AttributeError('run solve_f first')
	
	
	@property
	def f2a(self):
		
		if hasattr(self, 'F'):
			return round(self.F[1], 4)
			
		else:
			raise AttributeError('run solve_f first')
			
	
	@property
	def f2b(self):
		
		if hasattr(self, 'F'):
			return round(self.F[2], 4)
			
		else:
			raise AttributeError('run solve_f first')
	
	
	@property
	def f3(self):
		
		if hasattr(self, 'F'):
			return round(self.F[3], 4)
			
		else:
			raise AttributeError('run solve_f first')	
		

class CentCeq(BaseC):
	'''
	Central carbon atom with a singlet "s", a doublet with coupling constant "d",
	and a triplet "t"
	'''
	
	def __init__(self, name, I):
		'''
		Parameters
		name: str
			name of carbon atom, e.g. "Ala_alpha"
		I: array or list of float
			observed relative multiplet intensities of "s", "d", and "t"
		'''
		
		super().__init__(self, I)
		
		if self.I.size != 3:
			raise ValueError('The length of multiplet intensities of central carbon atom with equal coupling constants to the adjacent carbons should be 3')
	
	
	@property
	def f1(self):
		
		if hasattr(self, 'F'):
			return round(self.F[0], 4)
			
		else:
			raise AttributeError('run solve_f first')
	
	
	@property
	def f2(self):
		
		if hasattr(self, 'F'):
			return round(self.F[1], 4)
			
		else:
			raise AttributeError('run solve_f first')
			
	
	@property
	def f3(self):
		
		if hasattr(self, 'F'):
			return round(self.F[2], 4)
			
		else:
			raise AttributeError('run solve_f first')
	
	
	
	
