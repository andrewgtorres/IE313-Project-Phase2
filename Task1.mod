# This model is part of the IEMS 313 project in Spring 2021 at Northwestern
# University.
#
# AMPL model for the Security Constrained DC Optimal Power Flow problem.

###########################################################################
# Sets, defining the network structure
###########################################################################

  # Number of buses in the network
param num_buses;

  # Set of buses
set BUSES = 1..num_buses;

  # Set of branches (lines)
  #
  # The elements of the set are triplets (from bus, to bus, id)
  # The id is an integer that distinguishes different lines if there are
  # parallel ones.
set BRANCHES within (BUSES cross BUSES cross integer[1, Infinity));

  # Sets of generators
  #
  # The elements of the sets are pairs (bus, id)
  # The id is an integer that distinguishes different generators at the same bus.
set GENERATORS within (BUSES cross integer[1, Infinity));

###########################################################################
# Parameters of the network, specifically for the branches
###########################################################################

  # Reactances [in deg/MW]
param reactance{BRANCHES};
  # Line flow limits (from bus) [in MW]
param lower_line_limit{BRANCHES};
  # Line flow limits (to bus) [in MW]
param upper_line_limit{BRANCHES};
  # Susceptance (computed from reactances) [in MW/deg]
param b{(i,j,k) in BRANCHES} := 1/reactance[i,j,k];

###########################################################################
# Parameters for the generators
###########################################################################

  # Minimum power generation [in MW]
param pg_min{GENERATORS};
  # Maximum power generation [in MW]
param pg_max{GENERATORS};
  # Coefficient of linear cost function [in $/MWh]
param pg_cost{GENERATORS};
  # Participation factor [dimensionless]
param alpha{GENERATORS};

###########################################################################
# Parameters defining a demand scenario
###########################################################################

  # Demand at each bus [in MW]
  # If none is explicitly specified, the demand is zero.
param demand{BUSES}, default 0;

###########################################################################
# Specifications of contingencies
###########################################################################

  # Set of branch contingencies
set BRANCH_CONTINGENCIES within BRANCHES, default {};

  # Set of generator contingencies
set GENERATOR_CONTINGENCIES within GENERATORS, default {};

###########################################################################
# Optimization variables
###########################################################################

  # Power generation in nominal case [in MW]
var pg{(i,k) in GENERATORS}, >= pg_min[i,k], <= pg_max[i,k];
  # Voltage angles in nominal case [in deg]
var delta{BUSES};

  # Power generation in branch contingencies [in MW]
var pg_bc{(i,k) in GENERATORS, (ic, jc, kc) in BRANCH_CONTINGENCIES},
     >= pg_min[i,k], <= pg_max[i,k];
  # Voltage angles in branch contingencies
var delta_bc{i in BUSES, (ic, jc, kc) in BRANCH_CONTINGENCIES};

  # Power generation in generator contingencies [in MW]
var pg_gc{(i,k) in GENERATORS, (ic, kc) in GENERATOR_CONTINGENCIES: i != ic or k != kc},
     >= pg_min[i,k], <= pg_max[i,k];
  # Voltage angles in generator contingencies
var delta_gc{i in BUSES, (ic, kc) in GENERATOR_CONTINGENCIES};
  # Lost generation in generator contingency [in MW]
var Omega_gc{(ic, kc) in GENERATOR_CONTINGENCIES};

###########################################################################
# Objective function
###########################################################################

  # Minimize the total generation cost in nominal case
minimize total_cost: sum{(i,k) in GENERATORS} pg_cost[i,k]*pg[i,k];

###########################################################################
# Power balance equations
###########################################################################

  # Power balance in nominal case
subject to power_balance{i in BUSES}:
    sum{(i,k) in GENERATORS} pg[i,k] - demand[i]
     =
    sum{(i,j,k) in BRANCHES} b[i,j,k] * (delta[i]-delta[j]) +
    sum{(j,i,k) in BRANCHES} b[j,i,k] * (delta[i]-delta[j]);

  # Power balance in branch congingency
subject to power_balance_bc{i in BUSES, (ic, jc, kc) in BRANCH_CONTINGENCIES}:
    sum{(i,k) in GENERATORS} pg_bc[i,k,ic,jc,kc] - demand[i]
     =
    sum{(i,j,k) in BRANCHES: i != ic or j != jc or k != kc}
      b[i,j,k] * (delta_bc[i,ic,jc,kc]-delta_bc[j,ic,jc,kc]) +
    sum{(j,i,k) in BRANCHES: i != ic or j != jc or k != kc}
      b[j,i,k] * (delta_bc[i,ic,jc,kc]-delta_bc[j,ic,jc,kc]);

  # Power generator in branch congingency
subject to power_balance_gc{i in BUSES, (ic, kc) in GENERATOR_CONTINGENCIES}:
    sum{(i,k) in GENERATORS: i != ic or k != kc} pg_gc[i,k,ic,kc] - demand[i]
     =
    sum{(i,j,k) in BRANCHES}
      b[i,j,k] * (delta_gc[i,ic,kc]-delta_gc[j,ic,kc]) +
    sum{(j,i,k) in BRANCHES}
      b[j,i,k] * (delta_gc[i,ic,kc]-delta_gc[j,ic,kc]);

###########################################################################
# Line flow limits
###########################################################################

  # Line flow limit in nominal case
subject to line_flow_limit{(i,j,k) in BRANCHES}:
       lower_line_limit[i,j,k]
    <= b[i,j,k] * (delta[i]-delta[j])
    <= upper_line_limit[i,j,k];

  # Line flow limit in branch contingency
subject to line_flow_limit_bc{(i,j,k) in BRANCHES,
                              (ic, jc, kc) in BRANCH_CONTINGENCIES}:
       lower_line_limit[i,j,k]
    <= b[i,j,k] * (delta_bc[i,ic,jc,kc]-delta_bc[j,ic,jc,kc])
    <= upper_line_limit[i,j,k];

  # Line flow limit in generator contingency
subject to line_flow_limit_gc{(i,j,k) in BRANCHES,
                              (ic, kc) in GENERATOR_CONTINGENCIES}:
       lower_line_limit[i,j,k]
    <= b[i,j,k] * (delta_gc[i,ic,kc]-delta_gc[j,ic,kc])
    <= upper_line_limit[i,j,k];

###########################################################################
# Generator response in contingency
###########################################################################

  # Generator response in branch contingency.
  # All generators produce the same as in base case, no generation was lost
subject to response_bc{(i,k) in GENERATORS,
                       (ic, jc, kc) in BRANCH_CONTINGENCIES}:
    pg_bc[i,k,ic,jc,kc] = pg[i,k];

  # Generator response in generator contingency.
subject to response_gc{(i,k) in GENERATORS,
                       (ic, kc) in GENERATOR_CONTINGENCIES: i != ic or k != kc}:
    pg_gc[i,k,ic,kc] = pg[i,k] + alpha[i,k] * Omega_gc[ic,kc];
