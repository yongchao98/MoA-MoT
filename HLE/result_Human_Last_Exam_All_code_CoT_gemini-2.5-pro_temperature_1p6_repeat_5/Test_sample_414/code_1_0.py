import numpy as np
from scipy.optimize import fsolve

def solve_radiochemistry_problem():
    """
    Calculates the time between separation and first analysis for a Ba-140/La-140 sample.
    """
    # Half-lives in days
    T_half_Ba = 12.75
    T_half_La = 40.28 / 24.0  # Convert hours to days

    # Decay constants (lambda) in days^-1
    lambda_Ba = np.log(2) / T_half_Ba
    lambda_La = np.log(2) / T_half_La

    # The problem provides two activity measurements, 14 days apart.
    # A1 = 1.4 kBq/mL
    # A2 = 2.1 kBq/mL
    # The ratio of activities is A2 / A1
    activity_ratio = 2.1 / 1.4

    # We assume the measured activity is dominated by the daughter, La-140.
    # The activity of the daughter nuclide (La-140) at time t after separation is:
    # A_La(t) = C * (exp(-lambda_Ba * t) - exp(-lambda_La * t))
    # We need to solve for t1 in the equation: A_La(t1 + 14) / A_La(t1) = 1.5

    def equation_to_solve(t1):
        """
        Represents the equation: A_La(t1 + 14) / A_La(t1) - 1.5 = 0
        """
        # Numerator of the ratio: A_La(t1 + 14)
        num = np.exp(-lambda_Ba * (t1 + 14)) - np.exp(-lambda_La * (t1 + 14))

        # Denominator of the ratio: A_La(t1)
        den = np.exp(-lambda_Ba * t1) - np.exp(-lambda_La * t1)

        # Avoid division by zero if t1 is very close to 0
        if np.abs(den) < 1e-12:
            return np.inf

        return num / den - activity_ratio

    # Use a numerical solver to find the root of the equation.
    # An initial guess of 1 day is reasonable.
    initial_guess = 1.0
    t1_solution_days = fsolve(equation_to_solve, initial_guess)[0]

    # Print the equation being solved with the values plugged in
    print("The problem is solved by finding 't' in the following equation, which represents the ratio of La-140 activities 14 days apart:\n")
    print(f"(exp(-{lambda_Ba:.4f} * (t + 14)) - exp(-{lambda_La:.4f} * (t + 14))) / (exp(-{lambda_Ba:.4f} * t) - exp(-{lambda_La:.4f} * t)) = {activity_ratio:.1f}\n")
    print("Solving for 't' gives the time between the chemical separation and the first analysis.")
    print(f"The calculated time is approximately {t1_solution_days:.2f} days.")
    print("\nThis represents the time elapsed between purifying the Barium-140 and performing the first measurement.")
    print(f"\nFinal Answer:")
    print(f"<<<{t1_solution_days:.2f}>>>")

solve_radiochemistry_problem()