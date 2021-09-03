#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


r'''
python C:\Users\cwu\Desktop\Software\NMRratios\nmrratios\test.py
'''


import sys
sys.path.append(r'C:\Users\cwu\Desktop\Software\NMRratios')
from nmrratios import NMRrun, TermC, CentCne, CentCeq, Ratio


if __name__ == '__main__':
	
	# instantiate NMRratios to calculate transformation matrix K
	run1 = NMRrun(P1 = 0.08)
		
	# terminal carbon
	Ala_beta = TermC('Ala_beta', [0.18, 0.82])   # intensities will be standardized automatically
	f = Ala_beta.solve_f(run1.K_termC, standardize = True)
	print('Ala_beta', f)
	
	# central carbon with unequal coupling constant
	Ala_alpha = CentCne('Ala_alpha', [0.17, 0.09, 0.01, 0.73])
	Ala_alpha.solve_f(run1.K_centC_ne)
	print('Ala_alpha', Ala_alpha.f1, Ala_alpha.f2a, Ala_alpha.f2b, Ala_alpha.f3)
	
	# central carbon with equal coupling constant
	Ile_gamma1 = CentCeq('Ile_gamma1', [0.49, 0.43, 0.08])
	f = Ile_gamma1.solve_f(run1.K_centC_eq)
	

	# instantiate another NMRratios to have different K
	run2 = NMRrun(P1 = 0.12)
		
	# central carbon with unequal coupling constant
	Met_alpha = CentCne('Met_alpha', [0.16, 0.20, 0.26, 0.38])
	Met_alpha.solve_f(run2.K_centC_ne)
	#print(Met_alpha.f1, Met_alpha.f2a, Met_alpha.f2b, Met_alpha.f3)
	
	
	# list available ratios
	print(Ratio.available_ratios())
	
	# use built-in ratio
	from nmrratios import PEP_from_OAA
	
	print(PEP_from_OAA)
	
	value = PEP_from_OAA.evaluate(Asp_alpha_f2b = 0.5, Phe_alpha_f2b = 0.1)
	print('PEP_from_OAA', value)
	
	# redefine the formula
	PEP_from_OAA.redefine('Tyr_alpha_f2b/Asp_alpha_f2b', ['Tyr_alpha_f2b', 'Asp_alpha_f2b'])
	
	print(PEP_from_OAA)
	
	value = PEP_from_OAA.evaluate(Asp_alpha_f2b = 0.5, Tyr_alpha_f2b = 0.2)
	print('new PEP_from_OAA', value)
	
	# self-defined ratios
	myRatio = Ratio('myRatio', 'x**2/(x+y)', ['x', 'y'])
	value = myRatio.evaluate(y = 0.6, x = 0.2)
	print(myRatio)
	print('myRatio', value)
	
	
	
	
	
	
	
	
	
	
	
		
		
		