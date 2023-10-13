
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# try using on my data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
import pandas as pd
import gurobipy as gp
from gurobipy import GRB
import logging

# Create an environment with your WLS license
license_file = '/Users/alexgorelick/gurobi.lic'
with open(license_file) as file:
    lines = [line.rstrip() for line in file]

params = {
        "WLSACCESSID": lines[3].split('=')[1],
        "WLSSECRET": lines[4].split('=')[1],
        "LICENSEID": int(lines[5].split('=')[1]),
        }

# Read the data from 'mydata.txt' into a pandas DataFrame
#dat = pd.read_csv('/Users/alexgorelick/lab_repos/lpASCN/C161_data_for_gurobi_long.tsv', delimiter='\t')
dat = pd.read_csv('/Users/alexgorelick/Dropbox (Partners HealthCare)/MGH CLINICAL/ALEX/germline_data/lpASCN/MDA16/processed_data/MDA16_q=10_Q=20_P=0_1000kbs_data_for_gurobi_long_no_Liv1_LN1.tsv', delimiter='\t')
threshold = 0.02864416 ## based on data after removing LN1 and Liv2. 20th percentile of abs(logR(snp) - logR(segment))

# get lengths of variables
Samples = dat['sample'].unique()
Segments = dat['segment'].unique()
n_Samples = len(Samples)
dat.set_index(['sample','segment'], inplace=True) ## set 3 indices, sample, dipLogR, segment

env = gp.Env(params=params)
model = gp.Model(env=env)

m = {}
b = {}
for t in Samples:
    m[t] = model.addVar(vtype=GRB.CONTINUOUS, name='m_'+str(t), lb=0.25, ub=4)
    b[t] = model.addVar(vtype=GRB.CONTINUOUS, name='b_'+str(t), lb=-1.0, ub=1.0)

# Create auxilary continuous variables err s.t. err{t,s} = |y{t,s} - (x_{t,s}*m_{t}-b{t})|
# Create binary variables match s.t. match{s}==1 if sum_s(err{t,s} <= threshold) == num_samples
err = {}
match = {}

for t in Samples:
    for s in Segments:
        err[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='err_'+str(t)+','+str(s))
        match[t, s] = model.addVar(vtype=GRB.BINARY, name='match_'+str(t)+','+str(s))
        silence = model.addConstr(err[(t, s)] >= dat.loc[t,s].y - (dat.loc[t,s].x*m[t]-b[t]), name='c1_'+str(t)+','+str(s))
        silence = model.addConstr(err[(t, s)] >= -dat.loc[t,s].y + (dat.loc[t,s].x*m[t]-b[t]), name='c2_'+str(t)+','+str(s))
        model.addGenConstrIndicator(match[(t, s)], 1, err[(t, s)], GRB.LESS_EQUAL, threshold, name='c3_'+str(t)+','+str(s))


# create a variable for the objective function, which is a binary variable for whether a segment was a 'match' in at least 80% of samples.
allmatch = {}
for s in Segments:
    allmatch[s] = model.addVar(vtype=GRB.BINARY, name='allmatch_'+str(s))
    model.addGenConstrIndicator(allmatch[s], 1, gp.quicksum(match[(t, s)] for t in Samples), GRB.GREATER_EQUAL, 0.8*n_Samples, name='c4_'+str(s))


model.update()
model.setObjective(gp.quicksum(allmatch[s] for s in Segments), GRB.MAXIMIZE)
model.update()
model.optimize()


# Check if a solution is found
if model.status == GRB.OPTIMAL:
    print("Optimal Solutions Found: "+str(model.SolCount))
    df = pd.DataFrame([])
    for w in range(model.SolCount):
        print("Extracting values for optimal solution "+str(w))
        model.setParam(GRB.Param.SolutionNumber, w)

        # Create a dictionary to store the selected profiles for each tumor sample
        sample_m = {t: None for t in Samples}
        sample_b = {t: None for t in Samples}
        segment_allmatch = {s: None for s in Segments}

        for t in Samples:         
            m_varname='m_'+str(t)
            b_varname='b_'+str(t)
            m_value = model.getVarByName(m_varname).X
            b_value = model.getVarByName(b_varname).X
            sample_m[t] = m_value
            sample_b[t] = b_value

        for s in Segments:
            allmatch_varname='allmatch_'+str(s)
            allmatch_value = model.getVarByName(allmatch_varname).X
            segment_allmatch["same_"+str(s)] = allmatch_value

        out_m = pd.DataFrame.from_dict(sample_m, orient='index')
        out_m = out_m.transpose()
        out_m = out_m.set_axis("m_"+Samples, axis=1)
        out_b = pd.DataFrame.from_dict(sample_b, orient='index')
        out_b = out_b.transpose()
        out_b = out_b.set_axis("b_"+Samples, axis=1)
        out_am = pd.DataFrame.from_dict(segment_allmatch, orient='index')
        out_am = out_am.transpose()

        result = pd.concat([out_m, out_b, out_am], axis=1)
        result['w'] = w
        df = pd.concat([df, result], axis=0)

    df.to_csv('mb_values.tsv')

else:
    print("No optimal solution found.")


# Dispose of the Gurobi model
model.dispose()




