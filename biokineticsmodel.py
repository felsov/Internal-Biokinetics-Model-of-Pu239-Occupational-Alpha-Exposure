############################################ACTIVITY DISTRIBUTION COMPUTATION
import math
t=24 #amount of time since first exposure in hours

A_perday = 7.29 * 10**6 #activity of Pu239 inhaled (Bq) per mining shift
shift = 8 #shift length in hours
offshift = 16 #time not at wokr in hours
weekhours = 168 #hours in week
weekend = 48 #hours in weekend
lambda_pu = 0.693 / (24110 * 365) #radioacitve decay constant of Pu-239  in d^-1

#Decay constants:
lambda_LtB = 0.1000989 #clearance from lung to blood in d^-1
lambda_LtGI = 0.2   # clearance from distal bronchioles -> main bornchi (-> GI) d^-1
lambda_lung_total = lambda_pu + lambda_LtB + lambda_LtGI

lambda_LitB = 0.03 #clearance from liver to blood in d^-1
lambda_liver_total = lambda_LitB + lambda_pu

lambda_StB = 0.0076 #clearance from skeleton to blood in d^-1
lambda_skeleton_total = lambda_StB + lambda_pu

lambda_BtS = 0.3235 #clearance from blood to skeleton in d^-1
lambda_BtoLi = 0.1941 #clearance from blood to liver in d^-1
lambda_blood_total = lambda_StB + lambda_BtoLi + lambda_pu

stuck_fraction = 0.02  #the 2 percent of lung Pu that doesn't get cleared evidenced by literature
A_lung_stuck = 0
A_lung_clearing = 0
A_liver = 0
stuck_skeleton_fraction = 0.007  #the .7 percent of skeletal Pu that doesn't get cleared evidenced by literature
A_skeleton_stuck = 0
A_skeleton_clearing = 0
A_blood = 0

#For-loop to caluclate activity from differential equations:
dt = 1  #in hours
for hour in range(t):
    hour_of_the_week = hour % weekhours #take remainder to find what hour in the week it is at
    if hour_of_the_week >= (weekhours - weekend): #scenario of weekend (no work -> I(t) = 0)
        I_t = 0
    else:
        hour_of_the_day = hour_of_the_week% 24 #take remainder to find what hour it is at

        if hour_of_the_day < shift: #scenario where its either the hour in day of a shift or not
            I_t = A_perday / shift #average activity over the shift length
        else:
            I_t = 0

    A_lung_clearing = A_lung_clearing * math.exp(-lambda_lung_total * (dt / 24)) + I_t * (1 - stuck_fraction)
    A_lung_stuck = A_lung_stuck * math.exp(-lambda_pu * (dt / 24)) + I_t * stuck_fraction
    A_lung = A_lung_clearing + A_lung_stuck

    A_liver = A_liver * math.exp(-lambda_liver_total * (dt / 24)) + lambda_BtoLi * A_blood * (dt / 24)

    A_skeleton_clearing = A_skeleton_clearing * math.exp(-lambda_skeleton_total * (dt / 24)) + lambda_BtS * A_blood * (dt / 24) * (1 - stuck_skeleton_fraction)
    A_skeleton_stuck = A_skeleton_stuck * math.exp(-lambda_StB * (dt / 24)) + lambda_BtS * A_blood * (dt / 24) * stuck_skeleton_fraction
    A_skeleton = A_skeleton_clearing + A_skeleton_stuck

    A_blood = A_blood * math.exp(-lambda_blood_total * (dt / 24)) + lambda_LtB * A_lung * (dt / 24) + lambda_LitB * A_liver * (dt / 24) + lambda_StB * A_skeleton * (dt / 24)

print("Activity in lung, liver, skeleton, and blood in Bq:", A_lung, A_liver, A_skeleton, A_blood)

############################################ABSORBED DOSE COMPUTATION
time_passed = 28800 # seconds
alpha_E = 8.3897 * 10**(-13) #J/one alpha decay

lung_SAF = 4.8706
liver_SAF = 0.46615
skeleton_SAF = 1.19
blood_SAF = 0.001

D_lung = (A_lung * time_passed) * lung_SAF * alpha_E
D_liver = (A_liver * time_passed) * liver_SAF * alpha_E
D_skeleton = (A_skeleton * time_passed) * skeleton_SAF * alpha_E
D_blood = (A_blood * time_passed) * blood_SAF * alpha_E

print("Absorbed dose in lung, liver, skeleton, and blood in Gy:", D_lung, D_liver, D_skeleton, D_blood)

############################################EQUIVALENT DOSE COMPUTATION
eD_lung = D_lung * 20
eD_liver = D_liver *20
eD_skeleton = D_skeleton *20
eD_blood = D_blood*20
print("Equivlanet dose of lung, liver, skeleton, and blood in Sv:", eD_lung, eD_liver, eD_skeleton, eD_blood)

############################################EFFECTIVE DOSE COMPUTATION
E = (0.12)*eD_lung + (0.04)*eD_liver + (0)*eD_blood + (0.12)*eD_skeleton
print("The approximated effective dose in Sv:", E, "in time", t, "hours")
