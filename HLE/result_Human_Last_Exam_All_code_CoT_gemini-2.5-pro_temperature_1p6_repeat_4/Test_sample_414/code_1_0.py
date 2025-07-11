import numpy as np

def solve_radiochemistry_problem():
    """
    This script solves for the time between chemical separation and the first
    measurement based on the transient equilibrium of the Ba-140/La-140 pair.
    """
    # --- Define Constants ---
    # Half-lives in days
    T_HALF_BA = 12.75
    T_HALF_LA = 1.678

    # Time between the two measurements in days
    DELTA_T = 14.0

    # Measured activity ratio
    A1 = 1.4
    A2 = 2.1
    TARGET_RATIO = A2 / A1

    # --- Calculations ---
    # Calculate decay constants (lambda = ln(2) / T_half) in units of days^-1
    lambda_ba = np.log(2) / T_HALF_BA
    lambda_la = np.log(2) / T_HALF_LA

    # --- Define the function for the activity ratio ---
    def calculate_activity_ratio(t):
        """
        Calculates the ratio of La-140 activity at time (t + DELTA_T) to the activity at time t.
        This is derived from the Bateman equation for transient equilibrium.
        """
        # Numerator: exp(-lambda_Ba * (t + DELTA_T)) - exp(-lambda_La * (t + DELTA_T))
        numerator = np.exp(-lambda_ba * (t + DELTA_T)) - np.exp(-lambda_la * (t + DELTA_T))
        
        # Denominator: exp(-lambda_Ba * t) - exp(-lambda_La * t)
        denominator = np.exp(-lambda_ba * t) - np.exp(-lambda_la * t)
        
        # Avoid division by zero if t=0
        if denominator == 0:
            return np.inf
        
        return numerator / denominator

    # --- Print the problem setup ---
    print("The key radioactive pair is Ba-140 -> La-140 due to the short irradiation time.")
    print(f"Half-life of Ba-140: {T_HALF_BA} days (decay constant λ_Ba = {lambda_ba:.4f} days^-1)")
    print(f"Half-life of La-140: {T_HALF_LA} days (decay constant λ_La = {lambda_la:.4f} days^-1)")
    print("\nThe equation to solve for 't' (the time of the first measurement in days) is:")
    print(f"{TARGET_RATIO:.1f} = [exp(-{lambda_ba:.4f}*(t + {DELTA_T})) - exp(-{lambda_la:.4f}*(t + {DELTA_T}))] / [exp(-{lambda_ba:.4f}*t) - exp(-{lambda_la:.4f}*t)]")

    # --- Numerical Search for the Solution ---
    # We will search for the value of 't' that makes our function equal to the target ratio.
    best_t = None
    smallest_diff = np.inf

    # Search for t from 0.01 days to 10 days with a small step.
    for t_guess in np.arange(0.01, 10.0, 0.001):
        ratio = calculate_activity_ratio(t_guess)
        diff = abs(ratio - TARGET_RATIO)
        if diff < smallest_diff:
            smallest_diff = diff
            best_t = t_guess

    # --- Output the Final Answer ---
    if best_t is not None:
        t_in_hours = best_t * 24
        print("\nSolving the equation numerically yields:")
        print(f"The approximate time between separation and the first analysis is {best_t:.2f} days.")
        print(f"This is approximately {t_in_hours:.1f} hours.")
    else:
        print("\nCould not find a numerical solution in the specified range.")

solve_radiochemistry_problem()