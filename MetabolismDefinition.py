from jmetal.core.problem import FloatProblem
from jmetal.core.solution import FloatSolution

from cobra.flux_analysis import flux_variability_analysis
from cobra import Model, Reaction, Metabolite

import random

import csv
import logging

class MetabolismDefinition:
    pass

    def __init__(self, path_metabolites, path_reactions):
        
        self.model = Model('Dummy FBA_Problem Instance')                           
        self.metabolites = {}        
        self.create_reactions(path_metabolites, path_reactions)
                    
    def create_metabolites(self, path):
        with open(path, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                
                self.metabolites[row['variable_name']] = Metabolite(row['id'], name=row['name'], formula=row['formula'], compartment=row['compartment'])

        return self.metabolites
    
    def create_reactions(self, path_metabolites, path_reactions):
        self.create_metabolites(path_metabolites)   

        with open(path_reactions, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                reaction = Reaction(row['name'])
                reaction.name = row['name']
                reaction.lower_bound = float(row['lower_bound'])
                reaction.upper_bound = float(row['upper_bound'])

                metab = row['metabolites'].replace('[', '').replace(']','').split(",")
                coefi = row['coeficients'].replace('[', '').replace(']','').split(",")

                for m in range(len(metab)):
                    reaction.add_metabolites({ self.metabolites[metab[m]] : float(coefi[m]) })
                    
                self.model.add_reactions([reaction])
                
        return self.model

    def show_reactions(self):
        
        print("Reactions")
        print("---------")
        for x in self.model.reactions:
            print("%s : %s %f %f" % (x.id, x.reaction, x.lower_bound, x.upper_bound))

    def show_metabolites(self):
        print("")
        print("Metabolites")
        print("-----------")
        for x in self.model.metabolites:
            print('%9s : %s' % (x.id, x.formula))

    def change_bounds(self, reaction_id, lb, ub):
        reaction= self.model.reactions.get_by_id(reaction_id)
        reaction.upper_bound = 1000      
        reaction.lower_bound = -1000     
        reaction.upper_bound = ub
        reaction.lower_bound = lb

    def get_model(self):
        return self.model

    def get_metabolites(self):
        return self.metabolites
