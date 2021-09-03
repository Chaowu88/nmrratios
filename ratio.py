#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '08/29/2021'
__version__ = '1.0'




from sympy import Symbol, lambdify
from sympy.parsing.sympy_parser import parse_expr




class Ratio:
	
	def __init__(self, name, formula, variables):
		'''
		Parameters
		name: str
			ratio name
		formula: str
			expression of ratio
		variables: list of str
			str parsered as variables in expression
		'''
		self.name = name
		self.formula = formula
		self.variables = list(set(variables))
		self.varDict = {var: Symbol(var) for var in variables}
		
		self.expr = parse_expr(formula, local_dict = self.varDict)
		self.lamExpr = lambdify(self.variables, self.expr)
	
	
	def __repr__(self):
		
		return '%s: %s' % (self.name, self.expr)
	
	
	def evaluate(self, **var_values):
		'''
		Parameters
		var_values: dict
			values for variables in the ratio expression
		'''
		
		try:
			values = [var_values[var] for var in self.variables]
		except:
			raise ValueError('values shoule be provided for all variables in the ratio expression')
		
		self.exprValue = self.lamExpr(*values)
		
		return round(self.exprValue, 4)
	
	
	def redefine(self, formula, variables):
		'''
		Parameters
		formula: str
			new expression of ratio
		variables: list of str
			str parsered as variables in new expression
		'''
		
		self.formula = formula
		self.variables = list(set(variables))
		self.varDict = {var: Symbol(var) for var in variables}
		
		self.expr = parse_expr(formula, local_dict = self.varDict)
		self.lamExpr = lambdify(self.variables, self.expr)
		
		
	@classmethod
	def available_ratios(cls):
		
		cls.ratios = ['Ratio', 
					  'PEP_from_OAA', 
					  'OAA_from_PEP', 
					  'Pyr_from_Mal_lb', 
					  'Pyr_from_Mal_ub',
					  'Pyr_from_ED_lb',
					  'Pyr_from_ED_ub',
					  'OAA_exch_Fum']
		
		return cls.ratios		

	
	
	
PEP_from_OAA = Ratio('PEP_from_OAA', 
					 'Phe_alpha_f2b/Asp_alpha_f2b', 
					 ['Phe_alpha_f2b', 'Asp_alpha_f2b'])

OAA_from_PEP = Ratio('OAA_from_PEP', 
					 '(Asp_alpha_f2a + Asp_alpha_f3)/(Phe_alpha_f2a + Phe_alpha_f3)',
					 ['Asp_alpha_f2a', 'Asp_alpha_f3', 'Phe_alpha_f2a', 'Phe_alpha_f3'])
	
Pyr_from_Mal_lb = Ratio('Pyr_from_Mal_lb',
						'(Ala_alpha_f2b - Phe_alpha_f2b)/(1 - Phe_alpha_f2b)',
						['Ala_alpha_f2b', 'Phe_alpha_f2b'])

Pyr_from_Mal_ub = Ratio('Pyr_from_Mal_ub',
						'(Ala_alpha_f2b - Phe_alpha_f2b)/(Asp_alpha_f2b - Phe_alpha_f2b)',
						['Ala_alpha_f2b', 'Phe_alpha_f2b', 'Asp_alpha_f2b'])

Pyr_from_ED_lb = Ratio('Pyr_from_ED_lb',
					   '(Ala_alpha_f3 - Phe_alpha_f3 + Pyr_from_Mal_lb*(Phe_alpha_f3 - Asp_alpha_f3))/(G6P_C2_f3 - Phe_alpha_f3)',
					   ['Ala_alpha_f3', 'Phe_alpha_f3', 'Pyr_from_Mal_lb', 'Asp_alpha_f3', 'G6P_C2_f3'])
					   
Pyr_from_ED_ub = Ratio('Pyr_from_ED_lb',
					   '(Ala_alpha_f3 - Phe_alpha_f3 + Pyr_from_Mal_ub*(Phe_alpha_f3 - Asp_alpha_f3))/(G6P_C2_f3 - Phe_alpha_f3)',
					   ['Ala_alpha_f3', 'Phe_alpha_f3', 'Pyr_from_Mal_ub', 'Asp_alpha_f3', 'G6P_C2_f3'])					   

OAA_exch_Fum = Ratio('OAA_exch_Fum',
					 '2*Asp_beta_f3/(Asp_alpha_f3 + Asp_beta_f3)',
					 ['Asp_alpha_f3', 'Asp_beta_f3'])



		
