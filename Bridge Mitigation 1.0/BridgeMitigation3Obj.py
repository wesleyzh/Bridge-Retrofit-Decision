"""
Birdget Network Mitigation Resource Allocation, version 1.0

Author: Weili Zhang
Date: 10/11/2015

Input:
Nodes
Arcs

Output:
Pareto front of stochastic multi-objective optimization

"""


import pandas as pd
from gurobipy import *
import random
import sys
import math
from nsga2 import Solution
from nsga2 import NSGAII
import copy


def Stochastic_IPW(solution, source, target):
    
    """compute the mean number of IPWs between s and t"""
    
    global nodes, arcs, bridges, undirection_bridge, bridge_dic, num_scenarios,  hosptials, resiendents
    
    IPW_sum = 0
    fail_status = {}
    
    #Gurobi model set
    m = Model('Max_Network_FLow')
    m.setParam( 'OutputFlag', 0) 
    m.setParam( 'LogToConsole', 0 )
    m.setParam( 'LogFile', "" )   
    m.params.threads = 1
    m.params.NodefileStart = 0.5
    m.params.timeLimit = 60    
    
    arcs = tuplelist(arcs)
    
    #define variables
    flow = {}
    for i, j in arcs:
        flow[i,j] = m.addVar(name='flow_%s_%s' % (i, j), obj=0, lb = 0, ub=1)   #flow on arc i,j in infrastructure
        
    IndependentPathWays = m.addVar(name='IPW', obj=-1, lb=0)
    
    m.update()
    
    for k in range(num_scenarios):
        
        m.reset()
        for c in m.getConstrs():
            m.remove(c)
        m.update()
        
        t = 0
        for head, tail in undirection_bridge:
            random_num = random.random()
            if solution[t] == 1 and random_num < (bridge_dic['Probability'][head, tail, 'retrofit', 'E']):
                fail_status[head, tail] = 1
                fail_status[tail, head] = 1
            elif  solution[t] == 0 and random_num < (bridge_dic['Probability'][head, tail, 'pre-retrofit', 'E'] ):
                fail_status[head, tail] = 1
                fail_status[tail, head] = 1
            else:
                fail_status[head, tail] = 0
                fail_status[tail, head] = 0
            
            t += 1
    
        for i in nodes:
            if i == source:
                m.addConstr(quicksum(flow[a, b] for a, b in arcs.select('*', i))-
                            quicksum(flow[a, b] for a, b in arcs.select(i, '*'))
                            == IndependentPathWays )
            elif i == target:
                m.addConstr(quicksum(flow[a, b] for a, b in arcs.select('*', i))-
                            quicksum(flow[a, b] for a, b in arcs.select(i, '*'))
                            == -1*IndependentPathWays )
            else:
                m.addConstr(quicksum(flow[a, b] for a, b in arcs.select('*', i))-
                            quicksum(flow[a, b] for a, b in arcs.select(i, '*'))
                            == 0)                
        
        for head, tail in undirection_bridge:
            m.addConstr(flow[head, tail] <= (1-fail_status[head, tail]))
            m.addConstr(flow[tail, head] <= (1-fail_status[tail, head]))

        m.update()
        
        m.optimize()
        if  m.status == 2 or m.status == 9:
            IPW_k_s = -1 * m.objval
        
        else:
            print m.status
        
            
        IPW_sum += IPW_k_s
        
    return IPW_sum / float(num_scenarios)

class MitigationSolution(Solution):
    
    '''
    Solution for Bridge Network Mitigation.
    '''
    def __init__(self):
        '''
        Constructor.
        '''
        
        global undirection_bridge
        
        Solution.__init__(self, 3)
        
        #self.xmin = 0.0
        #self.xmax = 1.0
        
        for _ in range(len(undirection_bridge)):
            self.attributes.append(random.randint(0, 1))
            #self.attributes.append(0)
        
        #print self.attributes
        
        self.evaluate_solution()
        
    def evaluate_solution(self):
        '''
        Implementation of method evaluate_solution() for objective functions.
        '''
        global resiendents, hosptials, arc_dic, budget, undirection_bridge
        
        Total_Expected_EIPW = 0
        num_pair = 0
        
        resident_hospital_IPW = {}
        
        for resident in resiendents:
            for emergency in hosptials:
                num_pair += 1
                resident_hospital_IPW[resident, emergency] = Stochastic_IPW(self.attributes, resident, emergency)
                Total_Expected_EIPW += resident_hospital_IPW[resident, emergency]
                #print resident, emergency, resident_hospital_IPW[resident, emergency]
        
        total_cost = 0
        t = 0
        for i, j in undirection_bridge:
            total_cost += self.attributes[t] * arc_dic['Retrofit Cost'][i, j]
            t += 1
        
        self.objectives[0] = total_cost
                
        #budget_penalty = 0
        #if total_cost <= budget + 1:
            #budget_penalty = 0
        #else:
            #budget_penalty = 10000000 * (total_cost - budget)
        
        self.objectives[1] =  -1 * Total_Expected_EIPW / float(num_pair)  #first objective value, nsga2 by default is to minimize
        
        resident_pair = []
        for i in resiendents:
            for j in resiendents:
                if j != i and (i, j) not in resident_pair and (j, i) not in resident_pair:
                    resident_pair.append((i, j))
        
        difference = 0        
                    
        for i, j in resident_pair:
            EIPW_i = 0
            EIPW_j = 0
            for e in hosptials:
                EIPW_i += resident_hospital_IPW[i, e]
                EIPW_j += resident_hospital_IPW[j, e]
            
            difference += (EIPW_i - EIPW_j) ** 2
        
        self.objectives[2] = math.sqrt(difference) 
        
        print self.objectives[0], - 1 * self.objectives[1], self.objectives[2]
        
        
    def crossover(self, other):
        '''
        Crossover of MitigationSolution.
        '''
        global undirection_bridge
        
        child_solution = MitigationSolution()
        index = random.randint(0, len(undirection_bridge)-1)
        
        child_solution.attributes = self.attributes[:index] + other.attributes[index:]
        
        
        return child_solution
    
    def mutate(self):
        '''
        Mutation of MitigationSolution.
        '''
        global undirection_bridge
        
        index = random.randint(0, len(undirection_bridge)-1)
        
        self.attributes[index] = 1 - self.attributes[index]
        
        

seed = 100
random.seed(seed)

#read file as input ******************************************************************************
node_df = pd.DataFrame(pd.read_csv('node.csv', header='infer'))  #read node table
arc_df =  pd.DataFrame(pd.read_csv('arc.csv', header='infer'))  #read arc table
arc_df['Arc'] = zip(arc_df['From'], arc_df['To'])  #merge two columns as one tuple column as arc
bridge_df = pd.DataFrame(pd.read_csv('Bridges.csv', header='infer'))
bridge_df['Bridge_Status_Damage'] = zip(bridge_df['From'], bridge_df['To'], bridge_df['Status'], bridge_df['Damage State'])
bridge_df['Bridge'] = zip(bridge_df['From'], bridge_df['To'])


#pre-process the dataframe to some data types in python

#create dictionary
node_dic = node_df.set_index('NodeID').to_dict()
arc_dic = arc_df.set_index('Arc').to_dict()
bridge_dic = bridge_df.set_index('Bridge_Status_Damage').to_dict()

nodes = node_df['NodeID'].tolist()      #list of nodes
arcs = arc_df['Arc'].tolist()           #list of arcs
bridges = bridge_df['Bridge'].tolist()  #list of bridges

if (set(bridges) < set(arcs)) == False:
    print "Bridge indexes are out of the set of arcs"
    sys.exit(0)
    

hosptials = []     #list of hospitals
resiendents = []   #list of resiendents

for key, value in node_dic['Type'].items():
    if value == "hospital" or value == 'Hospital':
        hosptials.append(key)
    elif value == "resident" or value == "Resident":
        resiendents.append(key)
        


for head, tail in arc_dic['Bidirection'].keys():
    if (tail, head) not in arcs:
        arcs.append((tail, head))
        for key in arc_dic.keys():
            arc_dic[key][tail, head] =  arc_dic[key][head, tail]

undirection_bridge = copy.deepcopy(set(bridges))
undirection_bridge = list(undirection_bridge)

for head, tail in bridges:
    if (tail, head) not in bridges:
        bridges.append((tail, head))
            


#nsga2 parameters ***************************************
num_scenarios = 100

#compute the scope for budget
#max_budget = sum(arc_dic['Retrofit Cost'].values()) / 2.0
#budget = max_budget
ini_num_population = 50
search_num_population = 20
num_generation = 30

nsga2 = NSGAII(3, 0.1, 1.0)


P = []
for i in range(ini_num_population):
    P.append(MitigationSolution())

nsga2.run(P, search_num_population, num_generation)

csv_file = open('nsga2_out_3obj.csv', 'w')
csv_file.write("" + "Cost" + "," + "AverageEIPW" + "," + "DifferenceEIPW" + "\n")
    
for i in range(len(P)):
    csv_file.write("" + str(P[i].objectives[0]) + ", " + str(-1*P[i].objectives[1]) + ", " + str(P[i].objectives[2]) + "\n")
    
csv_file.close()