import math
from scipy.optimize import fsolve

def solve_radiochemistry_problem():
    """
    Solves the radiochemistry problem by modeling the decay of Ba-140 and ingrowth of La-140.
    """
    # Half-lives in days
    T_p = 12.75  # Ba-140 half-life
    T_d = 1.678   # La-140 half-life

    # Decay constants (lambda = ln(2)/T) in day^-1
    lambda_p = math.log(2) / T_p
    lambda_d = math.log(2) / T_d

    # Initial activity of pure Ba-140 parent at the time of first analysis (t=0) in kBq/mL
    A_p_initial = 1.4

    # Final total activity measured in kBq/mL
    A_total_final = 2.1

    # The equation for total activity A(t) at a time t after separating the pure parent is:
    # A(t) = A_p(t) + A_d(t)
    # A(t) = A_p(0)*exp(-lambda_p*t) + (lambda_d/(lambda_d-lambda_p))*A_p(0)*(exp(-lambda_p*t) - exp(-lambda_d*t))
    # We are given A(t) = 2.1 and A_p(0) = 1.4 and we need to solve for t.

    # Let's define the function to find the root of: f(t) = A_total(t) - 2.1 = 0
    def total_activity_equation(t):
        term1 = A_p_initial * math.exp(-lambda_p * t)
        ingrowth_factor = lambda_d / (lambda_d - lambda_p)
        term2 = ingrowth_factor * A_p_initial * (math.exp(-lambda_p * t) - math.exp(-lambda_d * t))
        return term1 + term2 - A_total_final

    # An initial guess for the time t. We know the activity increases, so it's relatively short.
    # The daughter activity peaks around 5.6 days, so the total activity peak is also nearby.
    initial_guess = 3.0
    
    # Use a numerical solver to find the time t
    time_to_reach_final_activity = fsolve(total_activity_equation, initial_guess)[0]

    # Print out the parameters and the final equation to be solved
    print("--- Parameters and Setup ---")
    print(f"Parent nuclide: Ba-140, Half-life (T_p) = {T_p:.2f} days")
    print(f"Daughter nuclide: La-140, Half-life (T_d) = {T_d:.3f} days")
    print(f"Parent decay constant (λ_p): {lambda_p:.4f} day⁻¹")
    print(f"Daughter decay constant (λ_d): {lambda_d:.4f} day⁻¹")
    print(f"Initial Parent Activity (A_p(0)) after separation: {A_p_initial} kBq/mL")
    print(f"Target Total Activity (A_total(t)): {A_total_final} kBq/mL\n")

    print("--- Solving the Equation for time 't' ---")
    print("The equation to solve is: A_total(t) = A_p(t) + A_d(t)")
    print("Where:")
    print("  A_p(t) = A_p(0) * exp(-λ_p * t)")
    print("  A_d(t) = [λ_d / (λ_d - λ_p)] * A_p(0) * [exp(-λ_p * t) - exp(-λ_d * t)]\n")

    print("We are solving the following equation for 't':")
    # To demonstrate the equation with numbers, we can calculate the pre-factor for the ingrowth term
    ingrowth_prefactor = lambda_d / (lambda_d - lambda_p)
    print(f"{A_total_final:.1f} = {A_p_initial:.1f} * exp(-{lambda_p:.4f} * t) + {ingrowth_prefactor:.4f} * {A_p_initial:.1f} * [exp(-{lambda_p:.4f} * t) - exp(-{lambda_d:.4f} * t)]\n")

    print("--- Result ---")
    print(f"The approximate time required for the total activity to reach {A_total_final} kBq/mL is {time_to_reach_final_activity:.2f} days.")
    print("This time corresponds to the time between the sample irradiation and the first analysis.")

solve_radiochemistry_problem()