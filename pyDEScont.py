import numpy as np
import scipy.stats as ss
import time
import multiprocessing as mp
from sysdef import System
from excelIO import select_file, read_source, write_output

class Simulation():
    """
      Class used to instantiate a simulation run.

      Key Attributes:
      DURATION    Duration of the simulated time, expressed in the same
                  Time Units (TU) of the flowrate.
                  DURATION shall be provided by the user.
                  
      TRANSITORY  Duration of the initial warm-up period, not considered for 
                  performance measurements (avg. flowrate, avg. buffer level).
                  Default: 10% of DURATION
                  
      clock       Simulation clock, updated after each event.
      
      num_iter    Counter of simulated events ("iterations").
                  Potentially, to be used as a termination condition instead
                  of "clock" or other purposes (e.g. progress monitoring).

      event_table List of events, used to determine the next event and update
                  the simulation clock.

      Other attributes are explained in specific functions.
    """
    def __init__(self, source):
        """
        Initializes simulation parameters.

        Parameters:
        source       Source file path, required to set DURATION.
        """
        self.DURATION = int(read_source(source, "SIMULATION",
                                        1, 1, rows_to_skip=2))
        self.TRANSITORY = int(self.DURATION/10)
        
        self.clock = 0
        self.prev_clock = 0
        self.num_iter = 0
        self.event_table = [] # scheduled event list

    def init_server(self, server, first_iter = True):
        """
        Initializes server attributes for the simulation (clock = 0).
        For the following events (first_iter == False), it resets current
        flowrate, in order to be correctly updated by set_current_flowrate.

        Operating Attributes:
        state          Current server state.
        
        curr_flowrate  Current flowrate, based the "state" variable.
        
        next_state     Next server state.
        
        operations     It is used to model operation-dependent failures,
                       i.e. if a server is slowed down, it will break down
                       later than if it worked at full speed.

        IN_EVENT_TABLE Flag if the server is already in the event_table.
                       It affects how the event is updated.

        Performance Attributes:
        STATE_TIME     Time spent by server in the current state over the
                       entire simulation.
                       
        STATE_FREQ     Ratio between STATE_TIME and overall simulated time. 

        TOT_FLOW       Total flow processed by server over the entire sim.

        AVG_FLOWRATE   TOT_FLOW divided by the overall simulated time.
        """
        # if True, initialize all system variables for the simulation
        if first_iter == True:
            server.state = 0
            server.curr_flowrate = float(server.flowrate[0,server.state])
            server.next_state = 0
            
            server.operations = 0
            
            server.IN_EVENT_TABLE = False

            server.STATE_TIME = np.zeros((server.num_states,1))
            server.STATE_FREQ = np.zeros((server.num_states,1))
            server.TOT_FLOW = 0
            server.AVG_FLOWRATE = 0
        # else, update only curr_flowrate
        else:
            s = server.state
            server.curr_flowrate = float(server.flowrate[0,s])

    def init_buffer(self, buffer):
        """
        Initializes buffer attributes for the simulation (clock = 0).

        Operating Attributes:
        level          Current buffer level. It can be considered the 
                       queue length for the following server.

        IN_EVENT_TABLE Flag if the buffer is already in the event_table.
                       It affects how the event is updated.

        Performance Attributes:
        AVG_LEVEL      Average buffer level during the overall simulated time.
        """
        buffer.level = 0
        buffer.IN_EVENT_TABLE = False
        buffer.AVG_LEVEL = 0 
        
    def go(self, system):
        """
        Starts the main simulation cycle, initializing system conditions
        at clock = 0. Then, it updates simulation clock until it reaches the
        required duration.

        Parameters:
        system     Object of System class initialized via sysdef module.
                   It includes input data about servers, buffers and overall
                   network topology.
        """
        self.max_states = 0

        # Initialization of system conditions at clock = 0
        for i in range(system.NUM_SERVERS):
            self.init_server(system.M[i])
            if system.M[i].num_states > self.max_states:
                self.max_states = system.M[i].num_states
            if (i < system.NUM_SERVERS-1):
                self.init_buffer(system.B[i])

        # Simulation cycle
        while self.clock < self.DURATION:
            cond1 = self.clock > self.TRANSITORY
            cond2 = self.prev_clock < self.TRANSITORY
            if cond1 == True and cond2 == True:
                self.TRANSITORY = self.prev_clock

            ### STEP 1: SET CURRENT FLOWRATES according to actual system
            ### conditions, i.e. if servers are blocked/starved 
            self.set_current_flowrate(system)

            ### STEP 2: UPDATE TABLE OF EVENTS
            ### STEP 2.A: UPDATE TABLE OF EVENTS - SERVER EVENTS
            ###
            ### Data recorded for each server event:
            ### Column 0: "M" ("Machine", equivalent to "server" event)
            ### Column 1: Server index
            ### Column 2: Current server state
            ### Column 3: Next server state
            ### Column 4: Next event (transition) time
            
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
            
            ### STEP 2: UPDATE TABLE OF EVENTS
            ### STEP 2.B: UPDATE TABLE OF EVENTS - BUFFER EVENTS
            ###
            ### Data recorded for each buffer event:
            ### Column 0: "B" ("Buffer")
            ### Column 1: Buffer index
            ### Column 2: Next buffer state ("Full", "Empty" or nothing)
            ### Column 3: Next event time
            for j in range(system.NUM_BUFFERS):
                if system.B[j].IN_EVENT_TABLE == False:
                    event = ["B", j+1, "", 0]
                    self.event_table.append(event)                
                    system.B[j].IN_EVENT_TABLE = True
                    
                self.update_events(j, "B", system)

            ### STEP 2.C: FIND EVENT WITH MINIMUM TIME (NEXT EVENT)
            next_event = min(self.event_table, key=lambda x: x[-1])
            self.event_table.remove(next_event)

            ### STEP 3: UPDATE SYSTEM STATE
            self.prev_clock = self.clock
            self.clock = next_event[-1] # move simulation clock forward

            ### STEP 3.A: CHANGE STATE OF OBJECTS != NEXT EVENT
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

            ### STEP 3.B: CHANGE STATE OF OBJECT AFFECTED BY EVENT
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
                
            ### last computations before next cycle....
            self.num_iter += 1
            if self.num_iter%10000==0:
                print("Run",self.ID,":",round(self.clock),end="\r")
            ### END OF MAIN LOOP

        self.calc_sim_results(system) # save simulation results
        print("Run", self.ID,":", self.num_iter,"iterations")

    def extrapolate_b(self, j, system):
        """
        Returns flowrates in and out a specific buffer.
        """
        u = system.B[j].upstream-1
        d = system.B[j].downstream-1
        
        flow_up = system.M[u].curr_flowrate * system.B[j].perc_up
        flow_down = system.M[d].curr_flowrate * system.B[j].perc_down

        return flow_up, flow_down

    def set_current_flowrate(self, system):
        """
        Resets current servers' flowrate according to buffers' states.
        """
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
        continue_check, r = True, 1
        while r<=2 or continue_check==True:
            # if check_blocking == false,
            # check if servers are "starved" by empty upstream buffers 
            check_blocking = False
            sequence = empty_buffers
            if r%2 == 0:
            # if true, check if servers are "blocked" by full downstream buffers
                check_blocking = True
                sequence = reversed(full_buffers)
            continue_check = self.check_slowdown(system, check_blocking, sequence)
            r += 1

    def check_slowdown(self, system, check_blocking, sequence):
        """
        Checks flowrate slowdown effects due to starving (check_blocking = False)
        or blocking (check_blocking = True) effects.
        """
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
            ## Cycle may need to be to be continued only if network
            ## has multiple branches.
            if (system.M[a].curr_flowrate < prev_flowrate
                and system.multi_branch == True):
                continue_cycle = True
        return continue_cycle

    def update_events(self, obj, obj_type, system, t=None):
        """
        Updates server/buffer events.
        """
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
        """
        Updates performance attributes of servers and buffers after each run.
        """
        ## evaluate only if warm-up time has ended
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
        """
        Copies system performance attributes from the System object to
        the Simulation itself.
        """
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

def parallel_run(r, sim_object, system, seed_array):
    """
    Assigns random seed to a Simulation object, then starts the
    specific simulation run.
    Returns a simulation object, including performance measures.

    Parameters:
    r              Index of simulation run.
    sim_object     Simulation object.
    system         System to be simulated.
    seed_array     Array of random number generation seeds.
    """
    np.random.seed(seed = seed_array[r])
    sim_object[r].go(system)

    return sim_object[r]

def overall_result(sim_data, output_file, multiple_branch):
    """
    Calculates and saves average system performance measures
    for all the simulation runs.

    Parameters:
    sim_data          Simulation object.
    output_file       Excel file where results are saved.
    multiple_branch   True/False. If FALSE (no multiple branches),
                      only one flowrate result is reported.
    """
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
    write_output(output_file,"SIM_RESULTS", results_array, 
                 index, columns, rows_to_skip=1)

def sim_statistics(sim_data, output_file, system, n_runs, runtime):
    """
    Calculates general statistics regarding the simulation itself.

    Parameters:
    sim_data          Simulation object.
    output_file       Excel file where results are saved.
    system            System to be simulated.
    n_runs            Nr. of independent simulation runs performed.
    runtime           Simulation duration, expressed in real time between
                      first and last run.
    """
    sim_time = 0
    transitory = 0
    for run in sim_data:
        sim_time += run.clock
        transitory += run.TRANSITORY
        
    ind = ["Number of runs:",
           "Simulated time per run (time units):",
           "Warm-up time per run (time units):",
           "Simulation runtime (sec):"]
    cols = ["Values"]
    data_array = np.empty((4,1))
    data_array[0,0] = n_runs
    data_array[1,0] = (sim_time-transitory)/n_runs
    data_array[2,0] = transitory/n_runs
    data_array[3,0] = runtime
    
    if system.multi_branch == True:
        k = system.NUM_SERVERS
    else:
        k = 1
    write_output(output_file,"SIM_RESULTS", data_array, ind, cols,
                 rows_to_skip=2+(k+system.NUM_BUFFERS)+1)
    
    
def main():
    n_cpu = mp.cpu_count()
    pool = mp.Pool(n_cpu)
    source = select_file() ## Opens "Select File" dialog box
    start = time.time()
    sys = System(source) ## system variables' initialization
    sim = []
    NUM_RUNS = int(read_source(source,"SIMULATION", 1, 1, rows_to_skip=1))
    np.random.seed(seed=1010)
    rand_seed = np.random.randint(1, 100000, NUM_RUNS)
    for r in range(NUM_RUNS):
        run = Simulation(source)
        sim.append(run)
        sim[r].ID = r+1
    # Starmap allows parallelization of multiple simulation runs
    sim = pool.starmap(parallel_run,
                       [(r, sim, sys, rand_seed) for r in range(NUM_RUNS)])
    pool.close()
    end = time.time()
    overall_result(sim, source, sys.multi_branch)
    print("Simulation time (sec):", end-start)
    sim_statistics(sim, source, sys, NUM_RUNS, end-start)

if __name__ == "__main__":
    main()
