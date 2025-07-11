import numpy as np

def solve_radiochemistry_problem():
    """
    Calculates the time between sample preparation and first analysis
    based on the decay kinetics of the Ba-140/La-140 pair.
    """

    # --- Given values ---
    # Half-life of Ba-140 in days
    t_half_ba = 12.75
    # Half-life of La-140 in days (40.28 hours / 24 hours/day)
    t_half_la = 1.678
    # Time between the two measurements in days
    delta_t = 14.0
    # First measured activity in kBq/mL
    A1 = 1.4
    # Second measured activity in kBq/mL
    A2 = 2.1

    # --- Calculations ---
    # Calculate decay constants (lambda = ln(2) / T_half)
    lambda_ba = np.log(2) / t_half_ba
    lambda_la = np.log(2) / t_half_la

    # Calculate the ratio of the two activities
    activity_ratio = A2 / A1

    # To find the time 't1', we solve the ratio equation:
    # R = [exp(-l_ba*(t1+dt)) - exp(-l_la*(t1+dt))] / [exp(-l_ba*t1) - exp(-l_la*t1)]
    # This can be rearranged to solve for t1 directly:
    # t1 = ln((R - exp(-l_la*dt)) / (R - exp(-l_ba*dt))) / (l_la - l_ba)
    
    exp_term_ba = np.exp(-lambda_ba * delta_t)
    exp_term_la = np.exp(-lambda_la * delta_t)
    
    numerator = np.log((activity_ratio - exp_term_la) / (activity_ratio - exp_term_ba))
    denominator = lambda_la - lambda_ba
    
    t1 = numerator / denominator

    # --- Output the result ---
    print("Step 1: Define physical constants.")
    print(f"Half-life of Ba-140: {t_half_ba} days")
    print(f"Half-life of La-140: {t_half_la} days")
    print(f"Decay constant for Ba-140 (λ_Ba): {lambda_ba:.4f} days^-1")
    print(f"Decay constant for La-140 (λ_La): {lambda_la:.4f} days^-1\n")
    
    print("Step 2: Set up the activity ratio equation from the two measurements.")
    print(f"First activity (A1): {A1} kBq/mL")
    print(f"Second activity (A2): {A2} kBq/mL")
    print(f"Time between measurements (Δt): {delta_t} days")
    print(f"Activity ratio (R = A2 / A1): {activity_ratio:.4f}\n")
    
    print("Step 3: Solve for the time 't' between separation and the first measurement.")
    print("The equation to solve is: R = [exp(-λ_Ba*(t+Δt)) - exp(-λ_La*(t+Δt))] / [exp(-λ_Ba*t) - exp(-λ_La*t)]")
    print(f"Solving this for 't' with the given values yields:")
    print(f"t = ln(({activity_ratio:.2f} - exp(-{lambda_la:.4f}*{delta_t:.0f})) / ({activity_ratio:.2f} - exp(-{lambda_ba:.4f}*{delta_t:.0f}))) / ({lambda_la:.4f} - {lambda_ba:.4f})")
    print(f"t ≈ {t1:.2f} days\n")

    print("The approximate time between the sample irradiation/separation and the first analysis is about 1 day.")
    
solve_radiochemistry_problem()
<<<1.03>>>