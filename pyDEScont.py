## DEScont performs a discrete-event simulation for systems where servers
## are modelled as continuous-time, discrete-state-space Markov chains with
## multiple up and down states. Feasible system layouts include not only
## lines, but even more complex topologies with multiple branches (useful for
## modelling assembly/disassembly manufacturing lines).
##
## THIS IS A BRAND NEW VERSION IN PYTHON 3 STARTED IN AUGUST 2020, USING
## CONCEPTS (BUT DIFFERENT APPROACH) TO A WORKING MATLAB VERSION MADE IN 2012.
##
## Main desired differences from the past:
## 1) MS Excel data initialization through pandas, instead of manual input
##    inside the script itself (or an obscure txt file);
## 2) Script optimization in order to save computation time;
## 3) Same data model (classes, modules and attributes) as the
##    analytical model, in order to use the same input files.
##
## For the time being, the material flow has only a single type.

import numpy as np
import scipy.stats as ss
import time
import multiprocessing as mp
from sysdef import System # SYSTEM INPUTS INITIALIZATION
from excelIO import select_file, read_source, write_output

##### N.B.: Servers and buffers' indexes in the source file are 1->N,
##### in line with scientific papers on these topic.
##### In Python, server M1 will be represented by index "0", up to N-1.
#####
##### In any case, source file index is stored in attribute "index",
##### but it won't be used in the code.

class Simulation():
    ### initialize single run simulation
    def __init__(self, source):
        self.DURATION = int(read_source("SIMULATION", 1, 1,
                                        source_file = source,
                                        rows_to_skip=1))
        self.TRANSITORY = int(self.DURATION/10)
        
        self.clock = 0
        self.prev_clock = 0
        self.num_iter = 0
        self.event_table = [] # scheduled event list

    def init_server(self, server, first_iter = True):
    # if first_iter = True, initialize all system variables for the simulation
        if first_iter == True:
            server.state = 0
            server.curr_flowrate = float(server.flowrate[0,server.state])
            server.next_state = 0
            
        # Used for Operation-dependent failures assumption
            server.operations = 0
            
            server.IN_EVENT_TABLE = False

            server.STATE_TIME = np.zeros((server.num_states,1))
            server.STATE_FREQ = np.zeros((server.num_states,1))
            server.TOT_FLOW = 0
            server.AVG_FLOWRATE = 0
    # else, update set only curr_flowrate
        else:
            s = server.state
            server.curr_flowrate = float(server.flowrate[0,s])

    def extrapolate_b(self, j, system):
            u = system.B[j].upstream-1
            d = system.B[j].downstream-1
            
            flow_up = system.M[u].curr_flowrate * system.B[j].perc_up
            flow_down = system.M[d].curr_flowrate * system.B[j].perc_down

            return flow_up, flow_down

    def check_slowdown(self, system, check_blocking, sequence):
        continue_cycle = False
        for j in sequence:
            flow_u, flow_d = self.extrapolate_b(j, system)
            BB_speed = min(flow_u, flow_d)
            if check_blocking == False:
                a = system.B[j].downstream-1
                percentage = system.B[j].perc_down
            else:
                a = system.B[j].upstream-1
                percentage = system.B[j].perc_up 
            prev_flowrate = system.M[a].curr_flowrate
            system.M[a].curr_flowrate = float(BB_speed/percentage)
            if system.M[a].curr_flowrate < prev_flowrate:
                continue_cycle = True
        return continue_cycle

    def set_current_flowrate(self, system):
        empty_buffers = []
        full_buffers = []
        for i in range(system.NUM_SERVERS):
            self.init_server(system.M[i], first_iter = False)

        for j in range(system.NUM_BUFFERS):
            if system.B[j].level == 0:
                empty_buffers.append(j)
            if system.B[j].level == system.B[j].capacity:
                full_buffers.append(j)

        ## in order to update current flowrate, alternatively check
        ## for starvation and blocking effects
        updated, r = True, 1
        while r<=2 or updated==True:
            # if check_blocking == false,
            # check if servers are "starved" by empty upstream buffers 
            check_blocking = False
            sequence = empty_buffers
            if r%2 == 0:
            # if true, check if servers are "blocked" by full downstream buffers
                check_blocking = True
                sequence = reversed(full_buffers)
            updated = self.check_slowdown(system, check_blocking, sequence)
            r += 1

    def update_events(self, obj, obj_type, system, t=None):

        obj_ID = obj+1

        for event in self.event_table:
            if (event[0] == "M" and obj_type == "M" and event[1] == obj_ID):
                event[4] = self.clock + t
            if (event[0] == "B" and obj_type == "B" and event[1] == obj_ID):
                flow_u, flow_d = self.extrapolate_b(obj_ID-1, system)
                diff = flow_u-flow_d

                if diff > 1e-6:
                    event[2] = "Full"
                    event[3] = float(self.clock + (system.B[obj_ID-1].capacity -
                                system.B[obj_ID-1].level)/(flow_u-flow_d))
                elif diff < -1e-6:
                    event[2] = "Empty"
                    event[3] = float(self.clock + system.B[obj_ID-1].level/
                                     (flow_d-flow_u))
                else:
                    event[2] = ""
                    event[3] = float("inf")

    def evaluate_statistics(self, obj, obj_type):
        if self.clock >= self.TRANSITORY:
            lti = self.clock - self.prev_clock # latest time interval
            TTI = self.clock - self.TRANSITORY # total time interval
            if obj_type == "M":
                s = obj.state
                obj.STATE_TIME[s] += lti
                obj.TOT_FLOW += obj.curr_flowrate * lti
                obj.AVG_FLOWRATE = obj.TOT_FLOW/TTI  
            if obj_type == "B":
                latest_avg_level = 0.5*(obj.level+obj.prec_level)*lti
                integral = obj.AVG_LEVEL*(TTI-lti)+latest_avg_level
                obj.AVG_LEVEL = integral/TTI
                
    def calc_sim_results(self, system):
        TTI = self.clock - self.TRANSITORY # total time interval
        self.avg_flowrate = np.zeros((system.NUM_SERVERS))
        self.state_freq = np.zeros((system.NUM_SERVERS, self.max_states))
        self.avg_level = np.zeros((system.NUM_BUFFERS))

        for i in range(system.NUM_SERVERS):
            self.avg_flowrate[i] = system.M[i].AVG_FLOWRATE
            for s in range(system.M[i].num_states):
                self.state_freq[i,s] = system.M[i].STATE_TIME[s]/TTI

        for j in range(system.NUM_BUFFERS):
            self.avg_level[j] = system.B[j].AVG_LEVEL
            
            
    ### initialize single run:
    ### for each server, state is set to the first one [0]
    ### for each buffer, level is set to zero.
    def go(self, system):
        self.max_states = 0
        for i in range(system.NUM_SERVERS):
            self.init_server(system.M[i])
            if system.M[i].num_states > self.max_states:
                self.max_states = system.M[i].num_states
            if (i < system.NUM_SERVERS-1):
                system.B[i].level = 0
                system.B[i].IN_EVENT_TABLE = False
                system.B[i].AVG_LEVEL = 0

        # Simulation cycle
        while self.clock < self.DURATION: #self.clock < self.DURATION:
            cond1 = self.clock > self.TRANSITORY
            cond2 = self.prev_clock < self.TRANSITORY

            if cond1 == True and cond2 == True:
                self.TRANSITORY = self.prev_clock

            self.set_current_flowrate(system)


            ### BLOCK B: SCHEDULED EVENT LIST GENERATION SECTION

            ### BLOCK B.1: UPDATE TABLE OF EVENTS - SERVER EVENTS
            ###
            ### GENERATE SERVER EVENT:
            ### IF EVENT ALREADY IN TABLE FOR THE SERVER, KEEP IT.
            ### WE ASSUME OPERATION-DEPENDENT FAILURES (ODF), I.E. IF MACHINE
            ### IS SLOWED DOWN, ITS EXPECTED TIME WILL BE POSTPONED.
            
            for i in range(system.NUM_SERVERS):

                if system.M[i].IN_EVENT_TABLE == False:
                    
                    s = int(system.M[i].state)

                    #find transition rates from current state
                    feasible_transition = np.copy(system.M[i].TM[s,:])

                    #trick to avoid "divide by zero" warning
                    feasible_transition[feasible_transition==0] = -1

                    # rework array to be workable by scipy.stats                    
                    feasible_transition = ((feasible_transition>0) *
                                           (1/feasible_transition))

                    feasible_transition = np.nan_to_num(feasible_transition,
                                                        copy=False)

                    # assumption: exponential time for all transitions
                    # N.B.: accurate for failures, less accurate for repairs
                    tt = ss.expon.rvs(scale=feasible_transition)

                    # time until next event, at nominal flowrate
                    time_min = np.min(tt[np.nonzero(tt)])
                    
                    # next server state according to earliest event 
                    next_s = int(np.where(tt == time_min)[0])
                    system.M[i].next_state = next_s

                    # operations: required time to next state if the server
                    # operates at nominal flowrate.
                    system.M[i].operations = time_min

                    # required time considering current flowrate
                    if system.M[i].flowrate[0,s] == 0:
                        duration = time_min
                    elif system.M[i].curr_flowrate == 0:
                        duration = float("inf")
                    else:
                        duration = float(system.M[i].operations *
                                     (system.M[i].flowrate[0,s] /
                                       system.M[i].curr_flowrate))                        

                    event = ["M", i+1, s, next_s, self.clock + duration]

                    self.event_table.append(event)
                    system.M[i].IN_EVENT_TABLE = True

                else:
                    # if already in event_table, update only its time
                    ## if server is down (flowrate = 0), don't change event time
                    s = int(system.M[i].state)
                    next_s = int(system.M[i].next_state)
                    
                    if system.M[i].flowrate[0,s] > 0:
                        if system.M[i].curr_flowrate == 0:
                            duration = float("inf")
                        else:
                            duration = float(system.M[i].operations *
                                     (system.M[i].flowrate[0,s] /
                                       system.M[i].curr_flowrate))

                        self.update_events(i, "M", system, duration)
            
            ### END OF BLOCK B.1

            ### END OF BLOCK B.2: UPDATE TABLE OF EVENTS - BUFFER EVENTS
            for j in range(system.NUM_BUFFERS):
                if system.B[j].IN_EVENT_TABLE == False:
                    event = ["B", j+1, "", 0]
                    self.event_table.append(event)                
                    system.B[j].IN_EVENT_TABLE = True
                    
                self.update_events(j, "B", system)

            ### END OF BLOCK B.2     

            # find event with minimum time
            next_event = min(self.event_table, key=lambda x: x[-1])
            self.event_table.remove(next_event)

            ### BLOCK C - SYSTEM STATE UPDATE
            self.prev_clock = self.clock
            self.clock = next_event[-1]

            ### BLOCK C.1 - CHANGE STATE OF OBJECTS != NEXT EVENT
            for event in self.event_table:
                a = event[1]-1
                if event[0] == "M":
                    s, next_s = event[2], event[3]
                    if system.M[a].flowrate[0,s] > 0:
                        nr_ops = float((self.clock - self.prev_clock) *
                                   (system.M[a].curr_flowrate /
                                     system.M[a].flowrate[0,s]))
                        system.M[a].operations -= nr_ops
                        self.evaluate_statistics(system.M[a], event[0]) 
                                                    
                else:
                    flow_u, flow_d = self.extrapolate_b(a, system)
                    system.B[a].prec_level = system.B[a].level 
                    system.B[a].level += float((self.clock - self.prev_clock) *
                                               (flow_u - flow_d))
                    self.evaluate_statistics(system.B[a], event[0])
                    
            ### END OF BLOCK C.1

            ### BLOCK C.2 - CHANGE STATE OF OBJECT AFFECTED BY EVENT
            a = next_event[1]-1
            if next_event[0] == "M":
                s, next_s = next_event[2], next_event[3]
                system.M[a].state = next_s
                system.M[a].operations = 0
                system.M[a].IN_EVENT_TABLE = False
                self.evaluate_statistics(system.M[a], next_event[0])

            else:
                system.B[a].IN_EVENT_TABLE = False
                if next_event[2] == "Full":
                    system.B[a].level = system.B[a].capacity
                else:
                    system.B[a].level = 0
                system.B[a].prec_level = system.B[a].level
                self.evaluate_statistics(system.B[a], next_event[0])

            ### END OF BLOCK C.2
            self.num_iter += 1
            if self.num_iter%10000==0:
                print("Run",self.ID,":",round(self.clock),end="\r")

            ### END OF BLOCK C
##            if self.clock == self.prev_clock:
##                #print("Warning! Simulation", self.ID, "stuck at time", self.clock)                    
##                for i in range(system.NUM_SERVERS):
##                    print("Vel.", i+1,":",system.M[i].curr_flowrate)
##                for j in range(system.NUM_BUFFERS):
##                    print("Livello", j+1,":",system.B[j].level)
        self.calc_sim_results(system)
        print("Run", self.ID,":", self.num_iter,"iterations")
        #print("Run", self.ID, "completed!")

def parallel_run(r, sim_object, system, seed_array):
    np.random.seed(seed = seed_array[r])
    sim_object[r].go(system)

    return sim_object[r]

def overall_result(sim_data, output_file, multiple_branch):
    ### results structure (example):
    ###  Values    | Average | CI 95% (±)
    ###   Flowrate | 0.822   | 0.021
    ###   Buffer 1 | 10.23   | 0.63
    ###   ...
    ###   Buffer N | 3.21    | 0.33
    
    n_runs = len(sim_data)
    n_servers = len(sim_data[0].avg_flowrate)
    n_buffers = len(sim_data[0].avg_level)
    avg_flowrate_array = np.zeros((n_servers,n_runs))
    avg_level_array = np.zeros((n_buffers,n_runs))

    row_nr = 0
    index = []

    for r in range(n_runs):
        avg_flowrate_array[:,r] = sim_data[r].avg_flowrate
        avg_level_array[:,r] = sim_data[r].avg_level

    if multiple_branch == False:
        results_array = np.zeros((1+n_buffers,2))
        index.append("Flowrate:")
        avg_flowrate_last = avg_flowrate_array[-1,:] 
        results_array[row_nr,0] = np.mean(avg_flowrate_last, axis=0)
        # Confidence interval of flowrate and buffer level mean values
        flowrate_CI95 = ss.bayes_mvs(avg_flowrate_last, alpha=0.95)[0][1]
        results_array[row_nr,1] = 0.5*(flowrate_CI95[1]-flowrate_CI95[0])
        row_nr += 1
    else:
        results_array = np.zeros((n_servers+n_buffers,2))
        for i in range(n_servers):
            string = "Server " + str(i+1) + " flowrate:"
            index.append(string)
            temp = ss.bayes_mvs(avg_flowrate_array[i,:], alpha=0.95)[0][1]
            results_array[row_nr,0] = np.mean(avg_flowrate_array[i,:])
            results_array[row_nr,1] = 0.5*(temp[1]-temp[0])
            row_nr += 1

    for j in range(n_buffers):
        string = "Buffer " + str(j+1) + " level:"
        index.append(string)
        temp = ss.bayes_mvs(avg_level_array[j,:], alpha=0.95)[0][1]
        results_array[row_nr,0] = np.mean(avg_level_array[j,:])
        results_array[row_nr,1] = 0.5*(temp[1]-temp[0])
        row_nr += 1
        
    columns = ["Average", "CI 95% (±)"]
    # write summary results to an excel file
    # Default: write to system input file
    write_output("SIM_RESULTS", results_array, index, columns,
                 output_file = output_file,
                 rows_to_skip=1)
    
def main():
    n_cpu = mp.cpu_count()
    pool = mp.Pool(n_cpu)
    source = select_file()
    start = time.time()
    sys = System(source) ## system variables' initialization
    sim = []
    NUM_RUNS = int(read_source("SIMULATION", 1, 1,source_file = source))
    np.random.seed(seed=1010)
    rand_seed = np.random.randint(1, 100000, NUM_RUNS)
    for r in range(NUM_RUNS):
        run = Simulation(source)
        sim.append(run) ## simulation initialization
        sim[r].ID = r+1
    sim = pool.starmap(parallel_run,
                       [(r, sim, sys, rand_seed) for r in range(NUM_RUNS)])
    pool.close()
    end = time.time()
    overall_result(sim, source, sys.multi_branch)
    print("Simulation time (sec):", end-start)
    # write simulation runtime statistics

    sim_time = 0
    transitory = 0
    for run in sim:
        sim_time += run.clock
        transitory += run.TRANSITORY
        
    ind = ["Number of runs:",
           "Simulated time per run (time units):",
           "Warm-up time per run (time units):",
           "Simulation runtime (sec):"]
    cols = ["Values"]
    data_array = np.empty((4,1))
    data_array[0,0] = NUM_RUNS
    data_array[1,0] = (sim_time-transitory)/NUM_RUNS
    data_array[2,0] = transitory/NUM_RUNS
    data_array[3,0] = end-start
    
    if sys.multi_branch == True:
        k = sys.NUM_SERVERS
    else:
        k = 1
    write_output("SIM_RESULTS", data_array, ind, cols,
                 output_file = source,
                 rows_to_skip=2+(k+sys.NUM_BUFFERS)+1)

if __name__ == "__main__":
    main()
