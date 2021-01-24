import numpy as np
from lib.excelIO import read_source, select_file

MAX_FLOWS = 100 ## Max number of inflows/outflows per server
MAX_STATES = 100 ## Max number of Markov states per server


class Server():
    def read_server_data(self, sheet, source):
        ## Server type: General (1 IN/1 OUT), Split (1 IN/N OUT),
        ## Merge (N IN/1 OUT)
        self.type = read_source(source, sheet, 1, 1, rows_to_skip=1)
        
        ## For Split/Merge servers, obtain percentage of flows for
        ## each downstream/upstream server (sum MUST BE 1)
        if(self.type == "Merge"):
            temp = read_source(source, sheet,
                               1, MAX_FLOWS, rows_to_skip=3)
            self.inflow = temp[np.isfinite(temp)]
            self.outflow = np.ones((1,1))
        elif(self.type == "Split"):
            temp = read_source(source, sheet,
                               1, MAX_STATES, rows_to_skip=3)
            self.inflow = np.ones((1,1))
            self.outflow = temp[np.isfinite(temp)]
        else:
            self.inflow = np.ones((1,1))
            self.outflow = np.ones((1,1))
            
        ## count nr. of inflows/outflows.
        ## Used later for assigning % of flow to the proper buffer. 
        self.num_inflows = len(self.inflow)
        self.num_outflows = len(self.outflow)
            
        ## Nominal flow rate for each state
        self.flowrate = read_source(source, sheet,
                                    1, MAX_STATES, rows_to_skip=5)
        self.flowrate = self.flowrate[self.flowrate>=0]

        ## Length of flowrate array is the number of server states
        self.num_states = len(self.flowrate)

        ## Transition Matrix of the Markov Chain representing the server
        self.TM = read_source(source, sheet,
                              self.num_states, self.num_states, rows_to_skip=8)

class Buffer():
    def __init__(self):
            self.index = None
            self.capacity = None
            self.upstream = None
            self.downstream = None

class System():
    def __init__(self,source):
        ## read all server data
        self.NUM_SERVERS = int(read_source(source,"NETWORK", 1, 1, rows_to_skip=1))
        self.M = []
        self.multi_branch = False
        for i in range(self.NUM_SERVERS):
            temp = Server()
            sheet_name = "M" + str(i+1)
            temp.read_server_data(sheet = sheet_name, source = source)
            self.M.append(temp)
            self.M[i].index = i+1
            if self.M[i].num_inflows > 1 or self.M[i].num_outflows > 1:
                self.multi_branch = True
        ## read all buffer data
        self.NUM_BUFFERS = self.NUM_SERVERS-1
        self.B = []
        B_mat = read_source(source,"NETWORK", 
                            self.NUM_BUFFERS, 4,
                            rows_to_skip=5)
        for j in range(self.NUM_BUFFERS):
            temp_B = Buffer()
            temp_B.index = j+1
            temp_B.capacity = B_mat[j][0]
            temp_B.upstream = B_mat[j][1]
            temp_B.downstream = B_mat[j][2]
            temp_B.pred, temp_B.succ = [], []
            self.B.append(temp_B)

        ## assign percentage of inflow/outflow to each buffer's upstream/
        ## downstream machine. Percentage assignment is based on buffers' order
        ## in the "NETWORK" input sheet.
        count_outflows = np.zeros((self.NUM_SERVERS,1),dtype=int)
        count_inflows = np.zeros((self.NUM_SERVERS,1),dtype=int)
        for j in range(self.NUM_BUFFERS):
            u = self.B[j].upstream-1 # reminder: M1 server has "0" index, etc.
            d = self.B[j].downstream-1            
            self.B[j].perc_up = self.M[u].outflow[count_outflows[u]]
            self.B[j].perc_down = self.M[d].inflow[count_inflows[d]]
            count_outflows[u] += 1
            count_inflows[d] += 1
            for k in range(self.NUM_BUFFERS):
                if(self.B[k].downstream == u+1):
                    self.B[j].pred.append(k)
                if(self.B[k].upstream == d+1):
                    self.B[j].succ.append(k)
            #print(j,":",self.B[j].pred, self.B[j].succ)
                    
        ## identify network topology (input/output ports)
        for i in range(self.NUM_SERVERS):
            count_in = np.count_nonzero(B_mat[:,1] == i+1)
            count_out = np.count_nonzero(B_mat[:,2] == i+1)
            
            if(count_in == 0):
                self.M[i].port = "OUT"
            elif(count_out == 0):
                self.M[i].port = "IN"
            else:
                self.M[i].port = ""

    # determine the nature of a buffer (input/output) according to the server
    # "port" attribute and determine buffers' update sequence for line
    # performance evaluation (either via simulation of analytical solution).
        self.SEQUENCE_arr = np.zeros(self.NUM_BUFFERS, dtype=int)
        num_inputs = 0
        num_outputs = 0
        count_inputs = 0
        count_outputs = 0
        count_others = 0        
        for j in range(self.NUM_BUFFERS):
            u = self.B[j].upstream-1 # reminder: M1 server has "0" index, etc.
            d = self.B[j].downstream-1         
            if(self.M[u].port == "IN"):
                self.B[j].port = "IN"
                num_inputs +=1
            elif(self.M[d].port == "OUT"):
                self.B[j].port = "OUT"
                num_outputs +=1               
            else:
                self.B[j].port = ""
        for j in range(self.NUM_BUFFERS):
            if(self.B[j].port == "IN"):
                self.SEQUENCE_arr[count_inputs] = j
                count_inputs += 1
            elif(self.B[j].port == "OUT"):
                self.SEQUENCE_arr[count_outputs-1] = j
                count_outputs -= 1
            else:
                self.SEQUENCE_arr[num_inputs+count_others] = j
                count_others += 1
        
        # convert SEQUENCE to list of integers for iterability
        self.SEQUENCE_arr.astype(int)
        self.SEQUENCE = self.SEQUENCE_arr.tolist()

if __name__ == "__main__":
    filename = select_file()
    sys = System(filename)
    for i in range(sys.NUM_SERVERS):
        print("Inflows, Server", i+1,":", sys.M[i].inflow)
        print("Outflows, Server", i+1,":", sys.M[i].outflow)
        print("Flowrate, Server", i+1,":", sys.M[i].flowrate)
        print("Transition Matrix, Server", i+1,":", sys.M[i].TM,"\n")
    for j in range(sys.NUM_BUFFERS):
        print("Upstream Server, Buffer", j+1,":", sys.B[j].upstream)
        print("Downstream Server, Buffer", j+1,":", sys.B[j].downstream)
        print("Capacity, Buffer", j+1,":", sys.B[j].capacity,"\n")


        

    
