
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# try using on my data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import pandas as pd
import gurobipy as gp
from gurobipy import GRB
import gurobipy_pandas as gppd

# Create an environment with your WLS license
license_file = '../gurobi.lic'
with open(license_file) as file:
    lines = [line.rstrip() for line in file]

params = {
        "WLSACCESSID": lines[3].split('=')[1],
        "WLSSECRET": lines[4].split('=')[1],
        "LICENSEID": int(lines[5].split('=')[1]),
        }

env = gp.Env(params=params)

# Create the model within the Gurobi environment
model = gp.Model(env=env)
#model.Params.PoolSearchMode = 2
#model.Params.PoolSolutions = 100

# Read the data from 'mydata.txt' into a pandas DataFrame
dat = pd.read_csv('processed_data/multifacets_autofit_cval=20_integerCN/dipLogR_sweep_aligned_segments_for_gurobi.tsv', delimiter='\t')

# get lengths of variables
Samples = dat['sample'].unique()
Profiles = dat['profile'].unique()
Segments = dat['segment'].unique()
dat.set_index(['sample','profile','segment'], inplace=True) ## set 3 indices, sample, dipLogR, segment

x = {}
for t in Samples:
    for p in Profiles:
        x[p, t] = model.addVar(vtype=GRB.BINARY, name='x_'+str(p)+','+str(t))

# Add constraint on x: Each tumor must be associated with exactly one profile
for t in Samples:
    silence = model.addConstr(gp.quicksum(x[(p, t)] for p in Profiles) == 1, name='c1_'+str(t))

# Create auxilary binary variables xbar{0,1} s.t. xbar_{p,t,pp,tt}=x_{p,t} * x_{pp,tt}.
xbar = {}
for t in Samples:
    print('sample=',t)
    for p in Profiles:
        for tt in Samples:
            for pp in Profiles:
                xbar[p, t, pp, tt] = model.addVar(vtype=GRB.BINARY, name='xbar_'+str(p)+','+str(t)+','+str(pp)+','+str(tt))
                silence = model.addConstr(xbar[(p, t, pp, tt)] <= x[p,t], name='c2_'+str(p)+','+str(t)+','+str(pp)+','+str(tt))
                silence = model.addConstr(xbar[(p, t, pp, tt)] <= x[pp,tt], name='c3_'+str(p)+','+str(t)+','+str(pp)+','+str(tt))
                silence = model.addConstr(xbar[(p, t, pp, tt)] >= x[p,t] + x[pp,tt] - 1, name='c4_'+str(p)+','+str(t)+','+str(pp)+','+str(tt))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create y[s] and add constraints on y
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create auxiliary binary variable z based on tmp1, tmp2, and create y{s} s.t. y_{s}=1 if s-th segment has the same value across all tumors and 0 otherwise
y = {}
for s in Segments:
    y[s] = model.addVar(vtype=GRB.BINARY, name='y_'+str(s))

for t in Samples:
    for p in Profiles:
        print('t:',t,',p:',p)
        for tt in Samples:
            for pp in Profiles:
                for s in Segments:
                    if int(dat.loc[t,p,s]) != int(dat.loc[tt,pp,s]) & int(dat.loc[tt,pp,s]) >= 0:
                        silence = model.addConstr(y[s] <= 1 - xbar[p,t,pp,tt], name='c5_'+str(p)+','+str(t)+','+str(pp)+','+str(tt)+','+str(s))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# do optimization
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

objective_expr = gp.quicksum(y[s] for s in Segments)

# Set objective
model.setObjective(objective_expr, GRB.MAXIMIZE)

# not sure what these lines do (are they necessary?)
model.update()  # process pending model modifications

# attempt to get an optimal solution
model.optimize()

# Print number of solutions stored
nSolutions = model.SolCount
print('Number of solutions found: ' + str(nSolutions))


# Check if an optimal solution has been found
if model.status == GRB.OPTIMAL:

    for w in range(model.SolCount):
        model.setParam(GRB.Param.SolutionNumber, w)

        # Create a dictionary to store the selected profiles for each tumor sample
        selected_profiles = {t: None for t in Samples}
        selected_segments = {s: None for s in Segments}

        for s in Segments:         
            segment_varname='y_'+str(s)    
            selected = model.getVarByName(segment_varname).X
            if selected > 0.5:
                selected_segments[s] = 1
            else:
                selected_segments[s] = 0

        # Extract the selected profiles from the solution
        for t in Samples:
            for p in Profiles:
                varname='x_'+str(p)+','+str(t)
                selected = model.getVarByName(varname).X
                if selected == 1:
                    selected_profiles[t] = p

        # Write the results to a tab-delimited file
        with open('selected_profiles_'+str(w)+'.tsv', 'w') as output_file:
            for t, p in selected_profiles.items():
                output_file.write(f"{t}\t{p}\n")

        # Write the results to a tab-delimited file
        with open('selected_segments_'+str(w)+'.tsv', 'w') as output_file:
            for s in selected_segments.items():
                output_file.write(f"{s}\n")

    print("Optimal solution found. Results written to selected_profiles.tsv.")
else:
    print("No optimal solution found.")


# Dispose of the model
model.dispose()

# compute IIS(), then export model.export('model')
