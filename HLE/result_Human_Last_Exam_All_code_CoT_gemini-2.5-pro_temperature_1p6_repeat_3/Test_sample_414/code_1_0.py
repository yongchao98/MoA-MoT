import numpy as np
from scipy.optimize import fsolve

def solve_radiochemistry_problem():
    """
    This script solves for the time between irradiation and analysis of a radioactive sample
    by modeling the parent-daughter decay of Ba-140 and La-140.
    """
    
    # Step 1: Define physical constants and measurements from the problem.
    # Half-lives in days
    T_half_Ba = 12.75
    T_half_La = 1.678
    
    # Measured activities in kBq/mL
    activity_1 = 1.4
    activity_2 = 2.1
    
    # Time between the two measurements in days
    delta_t = 14.0
    
    # Step 2: Calculate the decay constants (lambda = ln(2) / half-life).
    lambda_Ba = np.log(2) / T_half_Ba
    lambda_La = np.log(2) / T_half_La
    
    # The ratio of the two activities.
    activity_ratio = activity_2 / activity_1
    
    # Step 3: Define the function to be solved.
    # We need to find T such that the ratio of La-140 activities at T+delta_t and T equals activity_ratio.
    # A_La(t) is proportional to (exp(-lambda_Ba*t) - exp(-lambda_La*t)).
    # We are solving f(T) = 0 for the equation:
    # A(T+delta_t)/A(T) - activity_ratio = 0
    def equation_for_T(T):
        """
        Represents the equation we need to solve for T.
        T is the unknown time from irradiation to the first measurement.
        """
        # Activity is proportional to the term below
        activity_at_T = np.exp(-lambda_Ba * T) - np.exp(-lambda_La * T)
        activity_at_T_plus_14 = np.exp(-lambda_Ba * (T + delta_t)) - np.exp(-lambda_La * (T + delta_t))
        
        # Avoid division by zero if the initial guess is 0, though it's not a solution.
        if activity_at_T == 0:
            return float('inf')
            
        return (activity_at_T_plus_14 / activity_at_T) - activity_ratio

    # Step 4: Use a numerical solver to find T.
    # We provide an initial guess of 1.0 day.
    initial_guess_T = 1.0
    solution_T = fsolve(equation_for_T, initial_guess_T)[0]

    # Step 5: Print the final output as requested.
    print("Based on the decay dynamics of the Ba-140 / La-140 pair, we can determine the time since irradiation.")
    print("\nThe final equation being solved for the unknown time 'T' is:")
    
    equation_str = (
        f"{activity_2} / {activity_1} = "
        f"[exp(-{lambda_Ba:.4f} * (T + {delta_t:.0f})) - exp(-{lambda_La:.4f} * (T + {delta_t:.0f}))] / "
        f"[exp(-{lambda_Ba:.4f} * T) - exp(-{lambda_La:.4f} * T)]"
    )
    print(equation_str)
    
    print("\nSolving this equation gives the result:")
    print(f"\nThe approximate time between sample irradiation and the first analysis is {solution_T:.2f} days.")

solve_radiochemistry_problem()
<<<1.03>>>