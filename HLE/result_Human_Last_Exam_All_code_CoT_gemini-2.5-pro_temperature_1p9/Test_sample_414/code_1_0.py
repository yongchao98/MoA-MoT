import numpy as np
from scipy.optimize import fsolve

def solve_radiochemistry_problem():
    """
    Solves for the time between reactor irradiation and the first analysis
    of the radioactive sample.
    """

    # --- Step 1: Define Constants ---
    # Half-lives in days
    T_half_Ba = 12.75
    T_half_La = 1.678
    T_half_Sr = 28.8 * 365.25 # in days
    T_half_Y = 2.671

    # Decay constants (lambda = ln(2)/T_half) in days^-1
    lambda_Ba = np.log(2) / T_half_Ba
    lambda_La = np.log(2) / T_half_La
    lambda_Sr = np.log(2) / T_half_Sr
    lambda_Y = np.log(2) / T_half_Y

    # --- Step 2: Define activity growth functions ---
    # Ingrowth of La-140 from Ba-140 parent
    def f_Ba(t):
        if t <= 0: return 0
        return np.exp(-lambda_Ba * t) - np.exp(-lambda_La * t)

    # Ingrowth of Y-90 from Sr-90 parent (secular equilibrium approximation)
    def f_Y(t):
        if t <= 0: return 0
        # The Sr-90 decay is negligible over the experiment time
        return 1 - np.exp(-lambda_Y * t)

    # --- Step 3: Find t1 based on the key assumption ---
    # Assumption: A_La(t1) = A_Y(t1) = 0.7. The total activity change gives the constraint:
    # A_La(t1+14)/A_La(t1) + A_Y(t1+14)/A_Y(t1) = 3
    
    def equation_for_t1(t1):
        # Activity ratio for Ba/La system
        ratio_Ba = f_Ba(t1 + 14) / f_Ba(t1)
        # Activity ratio for Sr/Y system
        ratio_Y = f_Y(t1 + 14) / f_Y(t1)
        return ratio_Ba + ratio_Y - 3

    # Numerically solve for t1. An initial guess of 2 days is reasonable.
    t1_initial_guess = 2.0
    t1_solution = fsolve(equation_for_t1, t1_initial_guess)[0] # in days

    # --- Step 4: Calculate parent activities and T_cool ---
    # From A_La(t1)=0.7 and A_Y(t1)=0.7, find the parent activities at separation.
    # C_Ba is the activity normalization constant for the Ba/La system
    C_Ba_factor = lambda_La / (lambda_La - lambda_Ba)
    A_Ba_sep = 0.7 / (C_Ba_factor * f_Ba(t1_solution))

    # C_Sr is the activity of Sr-90 at separation
    A_Sr_sep = 0.7 / f_Y(t1_solution)
    
    # The ratio of parent activities depends on T_cool (cooling time)
    # A_Ba_sep / A_Sr_sep = R_initial * exp(-lambda_Ba * T_cool)
    
    # Calculate initial production ratio R_initial for short irradiation (e.g., 5 days)
    t_irr = 5 # days
    Y_Ba = 0.063 # Fission yield for Ba-140
    Y_Sr = 0.058 # Fission yield for Sr-90
    R_initial = (Y_Ba * (1 - np.exp(-lambda_Ba * t_irr))) / (Y_Sr * (1 - np.exp(-lambda_Sr * t_irr)))

    # Now solve for T_cool
    ratio_at_sep = A_Ba_sep / A_Sr_sep
    exp_term = ratio_at_sep / R_initial
    T_cool = -np.log(exp_term) / lambda_Ba

    # --- Step 5: Calculate the final answer ---
    total_time = T_cool + t1_solution

    print(f"Decay constant for Ba-140: {lambda_Ba:.4f} days⁻¹")
    print(f"Decay constant for La-140: {lambda_La:.4f} days⁻¹")
    print(f"Decay constant for Y-90:  {lambda_Y:.4f} days⁻¹")
    print("\nSolving for the time between separation and first measurement (t₁)...")
    print(f"Equation solved: Ratio_Ba(t₁) + Ratio_Sr(t₁) = 3")
    print(f"Result for t₁: {t1_solution:.2f} days")
    print("\nSolving for the cooling time before separation (T_cool)...")
    print(f"Initial production ratio A(Ba-140)/A(Sr-90): {R_initial:.1f}")
    print(f"Activity of Ba-140 at separation: {A_Ba_sep:.2f} kBq/mL")
    print(f"Activity of Sr-90 at separation: {A_Sr_sep:.2f} kBq/mL")
    print(f"Result for T_cool: {T_cool:.2f} days")
    print("\n--- Final Calculation ---")
    print("The approximate time between irradiation and the first analysis is T_cool + t₁.")
    print(f"Time = {T_cool:.2f} days + {t1_solution:.2f} days")
    print(f"Total Time = {total_time:.1f} days")

solve_radiochemistry_problem()