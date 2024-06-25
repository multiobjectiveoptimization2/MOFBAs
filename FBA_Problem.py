from pandas._libs.tslibs.offsets import delta_to_tick
from MetabolismDefinition import MetabolismDefinition

from jmetal.core.problem import FloatProblem
from jmetal.core.solution import FloatSolution

from cobra.flux_analysis import flux_variability_analysis
from cobra import Model, Reaction, Metabolite

import random

import csv
import logging

class FBA_Problem(FloatProblem):
    pass
        
    def __init__(self, vars, objs, cell: MetabolismDefinition= None, name = 'FBA_Problem'):
        super(FBA_Problem, self).__init__()
        
      
        self.model = cell.get_model()
        self.metabolites = cell.get_metabolites()
        
        self.solution = []
        self.reactions_lb = []
        self.reactions_ub = []
        
        for r in self.model.reactions:
            self.reactions_lb.append(r.lower_bound)
            self.reactions_ub.append(r.upper_bound)

        self.var_reactions = vars
        self.obj_reactions = objs

        self.number_of_variables = 2*len(vars)+2*len(objs)      
        self.number_of_objectives = len(objs)

        self.number_of_constraints = 1

        self.obj_directions = [self.MINIMIZE] * self.number_of_objectives 

        self.obj_labels = objs

        self.lower_bound = self.number_of_variables * [0.0]
        self.upper_bound = self.number_of_variables * [0.0]

        for i in range(self.number_of_variables):
            if i%2 == 0:
                if i < 2*len(vars):
                    x = self.model.reactions.get_by_id(vars[i//2])
                else:
                    x = self.model.reactions.get_by_id(objs[(i-2*len(vars))//2])
                lb = x.lower_bound      
                ub = x.upper_bound      
            else:
                lb = 0.0                                
                ub = 1.0

            self.lower_bound[i] = lb
            self.upper_bound[i] = ub

        
        
        logging.getLogger("cobra").setLevel(logging.ERROR)
        
    def compute_reference_set(self):
        r = []
        
        for obj_reaction in self.obj_reactions:
            s = self.external_fitness(obj_reaction)
            
            r.append(s)
        
        return r

    def change_bounds(self, reaction_id, lb, ub):
        reaction= self.model.reactions.get_by_id(reaction_id)      
        reaction.lower_bound = -1000
        reaction.upper_bound = 1000
        reaction.lower_bound = lb
        reaction.upper_bound = ub

    def external_fitness(self, r):
        self.model.objective = r    
        
        self.solution = self.model.optimize(r)    

        return self.solution        
        

    def fitness(self):                
        self.model.objective = self.obj_reactions[0]    #first in array leads the optimization
        
        self.solution = self.model.optimize()

        return self.solution        

    def evaluate(self, solution: FloatSolution) -> FloatSolution:
        
        x = solution.variables

        for i in range(len(self.var_reactions)+len(self.obj_reactions)):
            if i < len(self.var_reactions):
                r = self.model.reactions.get_by_id(self.var_reactions[i])                    
            else:
                r = self.model.reactions.get_by_id(self.obj_reactions[i-len(self.var_reactions)])
            
            r_idx = self.model.reactions.index(r.name)
            r_lb = self.reactions_lb[r_idx]
            r_ub = self.reactions_ub[r_idx]
            I = x[2*i]
            delta = x[2*i+1]

            

         
            lb = r_lb*delta
            
            if I < 0:
                ub = r_ub
            else:
                ub = I
                
            self.change_bounds(r.name, lb, ub)
           
                
        s = self.fitness()       

        if s.status == "optimal":
            for i in range(len(self.obj_reactions)):
                solution.objectives[i] = -1*self.solution.fluxes.get(self.obj_reactions[i])
        else:
            for i in range(len(self.obj_reactions)):
                solution.objectives[i] = 0            

        self.__evaluate_constraints(solution)

        return solution 

    def __evaluate_constraints(self, solution: FloatSolution) -> None:
        constraints = [0.0 for _ in range(self.number_of_constraints)]

        if self.solution.status != "optimal":
            constraints[0] = -1

        solution.constraints = constraints

        

    def get_name(self):
        return "FBA_MOOP"

