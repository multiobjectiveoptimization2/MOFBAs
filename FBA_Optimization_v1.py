from MetabolismDefinition import MetabolismDefinition

from jmetal.core.operator import R
from jmetal.core.problem import FloatProblem
from jmetal.core.solution import FloatSolution

from cobra.flux_analysis import flux_variability_analysis
from cobra import Model, Reaction, Metabolite

import random

import csv

import numpy as np
from pymoo.indicators.gd import GD
from pymoo.indicators.hv import HV
from pymoo.problems import get_problem

class FBA_Optimization_v1(FloatProblem):
    pass

    def __init__(self, obj_reactions, ref_point = [], cell: MetabolismDefinition= None, name = 'HV + DG + ObjReactions Optimization FBA Model'):
        super(FBA_Optimization_v1, self).__init__()
                      
        
        
        
       
        

       
        self.model = cell.get_model()
        self.metabolites = cell.get_metabolites()
        
        self.solution = []
        
        self.ref_point = np.array(ref_point)

        
        
        self.var_reactions = obj_reactions                  
                                                            
        self.obj_reactions = obj_reactions                  
        self.number_of_variables = len(obj_reactions)+1     
        self.number_of_objectives = 2 + len(obj_reactions)  

        self.number_of_constraints = 0                      

        obj_directions = []
        for i in range(0, self.number_of_objectives):
            obj_directions.append(self.MINIMIZE) 
                
        self.obj_directions = obj_directions

        self.obj_labels = obj_reactions
        
        
        Zr = self.compute_reference_set() 
        self.Zr = []
        for s in Zr:
            sol = []
            for r in self.var_reactions:
                flux = s.fluxes.get(r)
                sol.append(flux)
            
            self.Zr.append(sol)
            
        self.Zr = np.array(self.Zr)
        
        
        
        
        for i in range(0, len(self.var_reactions)):
            r = self.var_reactions[i]
            s = Zr[i]           
            
            reaction= self.model.reactions.get_by_id(r)
            
            flux = s.fluxes.get(r)
            
            reaction.lower_bound = -1*flux
            reaction.upper_bound = flux

        self.reactions_lb = []
        self.reactions_ub = []
        
        for r in self.model.reactions:
            self.reactions_lb.append(r.lower_bound)
            self.reactions_ub.append(r.upper_bound)

        
        self.lower_bound = self.number_of_variables * [0.0]
        self.upper_bound = self.number_of_variables * [0.0]

        for i in range(self.number_of_variables-1):
            x = self.model.reactions.get_by_id(obj_reactions[i])
            lb = x.lower_bound      
            ub = x.upper_bound      

            self.lower_bound[i] = lb
            self.upper_bound[i] = ub

        self.lower_bound[self.number_of_variables-1] = 0
        self.upper_bound[self.number_of_variables-1] = 1
        
      
        
        

        self.GD = GD(self.Zr)
        self.HV = HV(ref_point = self.ref_point)

    def compute_reference_set(self):
        r = []
        
        for obj_reaction in self.var_reactions:
            s = self.external_fitness(obj_reaction)
            
            r.append(s)
        
        return r
                
    def external_fitness(self, r):
        self.model.objective = r    
        
        self.solution = self.model.optimize(r)    

        return self.solution        
        
    def fitness(self, min_reaction):                
        self.model.objective = min_reaction    
            
        self.solution = self.model.optimize()

        return self.solution        

    def change_bounds(self, reaction_id, val):
        reaction= self.model.reactions.get_by_id(reaction_id)
        
        lb = -1*val
        ub = val
        
        if val < 0:
            lb = val
            ub = -1*val
            
        reaction.upper_bound = ub
        reaction.lower_bound = lb
        

    def evaluate(self, solution: FloatSolution) -> FloatSolution:
        
        
        
    
        nv = solution.number_of_variables-1
        delta = 1.0 / nv
        ini = 0
        fin = delta
        min_obj = -1
        val = solution.variables[nv]
        
        for i in range(0, nv):
            if val >= ini and val < fin:
                min_obj = i
                break
            ini = fin
            fin = fin + delta

        x = solution.variables       
        
        for i in range(0, len(self.var_reactions)):
            r = self.var_reactions[i]
            val = x[i]
            
            self.change_bounds(r, val)
            
        s = self.fitness(self.var_reactions[min_obj])       

        if s.status == "optimal":
            
            A = []
            for r in self.var_reactions:
                flux = s.fluxes.get(r)
                A.append(flux)
                
            A = np.array(A)
            H = self.HV(A)
            DG = self.GD(A)
            
            for i in range(0, self.number_of_objectives):
                if i == 0:
                    solution.objectives[i] = -1*H
                elif i == 1:
                    solution.objectives[i] = DG
                else:
                    r = self.var_reactions[i-2]
                    flux = s.fluxes.get(r)
                    solution.objectives[i] = -1*abs(flux)
        else:
            H = 0
            DG = 1000000
            
            for i in range(0, self.number_of_objectives):
                if i == 0:
                    solution.objectives[i] = -1*H
                elif i == 1:
                    solution.objectives[i] = DG
                else:
                    solution.objectives[i] = 0      
            
        return solution

    def __evaluate_constraints(self, solution: FloatSolution) -> None:
        return None

    def get_name(self):
        return "FBA_MOOP"



