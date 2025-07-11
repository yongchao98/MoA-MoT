import math

# --- Given Data ---
# Half-lives in days
T_half_Ba = 12.75
T_half_La = 1.678

# Measured activities in kBq/mL
A1 = 1.4
A2 = 2.1

# Time between measurements in days
dt = 14

# --- Calculations ---
# 1. Calculate decay constants (lambda) in units of days^-1
lambda_Ba = math.log(2) / T_half_Ba
lambda_La = math.log(2) / T_half_La

# 2. Calculate the ratio of the activities
activity_ratio = A2 / A1

# 3. Solve the equation for t1, the time of the first measurement.
# The equation is: t1 = ln( (ratio - exp(-lambda_Ba*dt)) / (ratio - exp(-lambda_La*dt)) ) / (lambda_Ba - lambda_La)
exp_term_Ba = math.exp(-lambda_Ba * dt)
exp_term_La = math.exp(-lambda_La * dt)

numerator_of_ln = activity_ratio - exp_term_Ba
denominator_of_ln = activity_ratio - exp_term_La

# Ensure the argument of the logarithm is positive
if numerator_of_ln <= 0 or denominator_of_ln <= 0 or (numerator_of_ln / denominator_of_ln) <=0:
    print("Error: Invalid arguments for logarithm. Cannot solve.")
else:
    log_argument = numerator_of_ln / denominator_of_ln
    t1 = math.log(log_argument) / (lambda_Ba - lambda_La)

    # --- Output Results ---
    print("This script solves for the time (t) of the first analysis using the relationship:")
    print("A2/A1 = [exp(-lambda_Ba*(t+dt)) - exp(-lambda_La*(t+dt))] / [exp(-lambda_Ba*t) - exp(-lambda_La*t)]")
    print("\nWhich rearranges to:")
    print("t = ln((A2/A1 - exp(-lambda_Ba*dt)) / (A2/A1 - exp(-lambda_La*dt))) / (lambda_Ba - lambda_La)")
    print("\n--- Input Values & Constants ---")
    print(f"Activity A1: {A1} kBq/mL")
    print(f"Activity A2: {A2} kBq/mL")
    print(f"Time between measurements (dt): {dt} days")
    print(f"Activity Ratio (A2/A1): {activity_ratio:.4f}")
    print(f"Half-life of Ba-140: {T_half_Ba} days")
    print(f"Half-life of La-140: {T_half_La} days")
    print(f"Decay constant lambda_Ba: {lambda_Ba:.6f} days^-1")
    print(f"Decay constant lambda_La: {lambda_La:.6f} days^-1")
    
    print("\n--- Equation Terms ---")
    print(f"A2/A1 - exp(-lambda_Ba*dt) = {activity_ratio:.4f} - {exp_term_Ba:.6f} = {numerator_of_ln:.6f}")
    print(f"A2/A1 - exp(-lambda_La*dt) = {activity_ratio:.4f} - {exp_term_La:.6f} = {denominator_of_ln:.6f}")
    print(f"Argument of ln: {log_argument:.6f}")
    print(f"lambda_Ba - lambda_La: {lambda_Ba - lambda_La:.6f} days^-1")
    
    print("\n--- Final Answer ---")
    print(f"The equation gives t = ln({log_argument:.6f}) / ({lambda_Ba - lambda_La:.6f})")
    print(f"The approximate time between irradiation and the first analysis is {t1:.2f} days.")

<<<1.03>>>