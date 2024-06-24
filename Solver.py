from FBA_Problem import FBA_Problem
from FBA_ProblemOptimized import FBA_ProblemOptimized
from FBA_Optimization_v1 import FBA_Optimization_v1
from FBA_Optimization_v2 import FBA_Optimization_v2
from MetabolismDefinition import MetabolismDefinition

from jmetal.algorithm.multiobjective import NSGAII
from jmetal.operator import SBXCrossover, PolynomialMutation

from jmetal.algorithm.multiobjective import MOEAD
from jmetal.operator import DifferentialEvolutionCrossover
from jmetal.util.aggregative_function import Tschebycheff

from jmetal.util.termination_criterion import StoppingByEvaluations
from jmetal.util.solution import get_non_dominated_solutions, print_function_values_to_file, print_variables_to_file

import csv

class Solver:
    def __init__(self):
        self.cell = None
        self.mo_fba = None
        self.reactions = None
        self.variables = None
        self.objectives = None
        self.algorithm = None
        self.results = None
        
    def get_optimization_elements(self, path_metabolites, path_reactions, obj_reactions):
        
        self.cell = MetabolismDefinition(path_metabolites, path_reactions)
    
        self.objectives = obj_reactions.split()
        self.reactions = []
        self.variables = []
        for r in self.cell.model.reactions:
            self.reactions.append(r.name)
        
            if r.name not in self.objectives:
                self.variables.append(r.name)    

        return self.reactions, self.variables, self.objectives

    def define_problem(self, problem_name):
        rp = [1000]*len(self.objectives)
        
        self.mo_fba = None
        if problem_name in 'FBA_Problem': 
            self.mo_fba = FBA_Problem(vars = self.variables, objs = self.objectives, cell=self.cell)        
        elif problem_name in 'FBA_ProblemOptimized': 
            self.mo_fba = FBA_ProblemOptimized(obj_reactions= self.objectives, ref_point = rp, cell=self.cell)
        elif problem_name in 'FBA_Optimization_v1': 
            self.mo_fba = FBA_Optimization_v1(obj_reactions= self.objectives, ref_point=rp, cell=self.cell)
        elif problem_name in 'FBA_Optimization_v2': 
            self.mo_fba = FBA_Optimization_v2(obj_reactions= self.objectives, cell=self.cell)
        else:
            print("Undefined Problem")
        
        return self.mo_fba
    
    def create_header(self):
        s = []

        s.append('configuration')
        s.append('experiment')
        s.append('name')
        s.append('status')

        for i in range(len(self.objectives)):
            s.append(self.objectives[i])

        

        for i in range(len(self.reactions)):
            s.append(self.reactions[i])

        return s

    def append_result(self, configuration, experiment, name, solution):
        s = []

        
        s.append(configuration)

        
        s.append(str(experiment))

        
        s.append(name)

        
        s.append(solution.status)

        
        for i in range(len(self.objectives)):
            s.append(solution.fluxes.get(self.objectives[i]))

        
        
        

        
        for i in range(1,len(self.reactions)+1):
            reaction_name = self.reactions[i-1]
            s.append(solution.fluxes.get(reaction_name));

        self.results.append(s)
    
    def define_algorithm(self, algorithm_name, population_size, max_evaluations):
        population_size = int(population_size)
        max_evaluations = int(max_evaluations)
        self.algorithm = None
        
        self.results = []

        
        rs = self.mo_fba.compute_reference_set()
        
        for i in range(0, len(rs)):
            s = rs[i]
            self.append_result(self.objectives[i], 0, self.objectives[i], s)
            
        if 'nsgaii' in algorithm_name:
            self.algorithm = NSGAII(
                problem=self.mo_fba,
                population_size=population_size,
                offspring_population_size=population_size,
                mutation=PolynomialMutation(probability=1.0 / self.mo_fba.number_of_variables, distribution_index=20),
                crossover=SBXCrossover(probability=1.0, distribution_index=20),
                termination_criterion=StoppingByEvaluations(max_evaluations=max_evaluations)
            )
        elif 'moead' in algorithm_name:
            neighbor_size = 10
            
            if population_size <= 20:
                neighbor_size = 4

            self.algorithm = MOEAD(
                problem=self.mo_fba,
                population_size=population_size,
                crossover=DifferentialEvolutionCrossover(CR=1.0, F=0.5, K=0.5),
                mutation=PolynomialMutation(probability=1.0 / self.mo_fba.number_of_variables, distribution_index=20),
                aggregative_function=Tschebycheff(dimension=self.mo_fba.number_of_objectives),
                neighbor_size=neighbor_size,
                neighbourhood_selection_probability=0.9,
                max_number_of_replaced_solutions=2,
                weight_files_path='_resources/',
                termination_criterion=StoppingByEvaluations(max_evaluations=max_evaluations)
            )
        elif 'fba' in algorithm_name:
            print('Done') 
        elif 'random' in algorithm_name: 
            self.algorithm = NSGAII(
                problem=self.mo_fba,
                population_size=population_size,
                offspring_population_size=population_size,
                mutation=PolynomialMutation(probability=1.0 / self.mo_fba.number_of_variables, distribution_index=20),
                crossover=SBXCrossover(probability=1.0, distribution_index=20),
                termination_criterion=StoppingByEvaluations(max_evaluations=population_size)
            )
        else:
            print('Undefined Algorithm')    

        return self.algorithm
    
    def solve(self, cfg):
        r, v, o = self.get_optimization_elements(cfg['#path_metabolites'], cfg['#path_reactions'],  cfg['#objective_reactions'])
        
        
        self.define_problem(cfg['#problem'])        
        
        
        num_exp = int(cfg['#runs'])
        
        if 'fba' not in cfg['#algorithm']:
            for run in range(1, num_exp+1):
                self.define_algorithm(cfg['#algorithm'], cfg['#population_size'], cfg['#max_evaluations'])
                self.algorithm.run()
                
                
                front = get_non_dominated_solutions(self.algorithm.get_result())

                
                id_name = cfg['#instance_name'] + '-' + cfg['#algorithm'] + '-' + cfg['#problem'] + '-' + cfg['#configuration_name'] + '-RUN_' + str(run)
                print_function_values_to_file(front, '_OUTPUT/FUN-' + id_name + '.txt')
                print_variables_to_file(front, '_OUTPUT/VAR-' + id_name + '.txt')

                
                for i in range(len(front)):
                    self.mo_fba.evaluate(front[i])
                    if self.mo_fba.solution.status == 'optimal':
                        self.append_result(cfg['#configuration_name'], run, 's'+str(i), self.mo_fba.solution)

    
                header = self.create_header()

                file_name = '_OUTPUT/RES-' + id_name + '.csv'

                with open(file_name, 'w', encoding='UTF8', newline='') as f:
                    writer = csv.writer(f)
        
                    writer.writerow(header)

                    writer.writerows(self.results)                
        else:
            self.define_algorithm(cfg['#algorithm'], cfg['#population_size'], cfg['#max_evaluations'])
          
            id_name = cfg['#instance_name'] + '-' + cfg['#algorithm'] + '-' + cfg['#problem'] + '-' + cfg['#configuration_name'] + '-RUN_1'
            
            header = self.create_header()

            file_name = '_OUTPUT/RES-' + id_name + '.csv'

            with open(file_name, 'w', encoding='UTF8', newline='') as f:
                writer = csv.writer(f)
        
                writer.writerow(header)

                writer.writerows(self.results)    
            
            f.close()
            
        return self.results
        
