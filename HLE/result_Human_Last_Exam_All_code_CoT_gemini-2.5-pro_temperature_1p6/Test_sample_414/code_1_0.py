import math

def solve_radiochemistry_problem():
    """
    Calculates the time between sample separation and the first analysis
    based on the decay of the Ba-140/La-140 pair.
    """
    # --- 1. Define constants ---
    # Half-lives in days
    T_Ba = 12.75  # days
    T_La_hours = 40.28 # hours
    T_La = T_La_hours / 24.0  # days

    # Decay constants (lambda = ln(2)/T) in days^-1
    lambda_Ba = math.log(2) / T_Ba
    lambda_La = math.log(2) / T_La

    # Activity measurements in kBq/mL
    A_La_1 = 1.4  # Activity at first measurement
    A_La_2 = 2.1  # Activity at second measurement

    # Time between measurements in days
    delta_t = 14.0

    # --- 2. Calculate intermediate values for the evolution equation ---
    # The evolution of the daughter (La-140) activity from measurement 1 to 2 is:
    # A_La_2 = A_La_1 * exp(-lambda_La*dt) + A_Ba_1 * [lambda_La/(lambda_La-lambda_Ba)] * [exp(-lambda_Ba*dt) - exp(-lambda_La*dt)]
    # We need to solve for A_Ba_1 (activity of Ba-140 at the first measurement).

    # Pre-calculate exponential terms
    exp_decay_La = math.exp(-lambda_La * delta_t)
    exp_decay_Ba = math.exp(-lambda_Ba * delta_t)

    # Pre-calculate the transient equilibrium constant K
    K = lambda_La / (lambda_La - lambda_Ba)

    # --- 3. Solve for the Ba-140 activity at the first measurement (A_Ba_1) ---
    # Rearrange the equation: A_Ba_1 = [A_La_2 - A_La_1*exp_decay_La] / [K * (exp_decay_Ba - exp_decay_La)]
    numerator = A_La_2 - A_La_1 * exp_decay_La
    denominator = K * (exp_decay_Ba - exp_decay_La)
    A_Ba_1 = numerator / denominator

    # --- 4. Solve for the time between separation and first analysis (t1) ---
    # At time t1 after separation, the ratio of activities is:
    # A_La_1 / A_Ba_1 = K * (1 - exp(-(lambda_La - lambda_Ba) * t1))
    # We can solve this for t1.

    # Rearrange to solve for the exponential term
    activity_ratio = A_La_1 / A_Ba_1
    exp_term = 1 - (activity_ratio / K)

    # Solve for t1
    # -(lambda_La - lambda_Ba) * t1 = log(exp_term)
    t1 = math.log(exp_term) / (-(lambda_La - lambda_Ba))
    
    # --- 5. Print the results step-by-step ---
    print("--- Analysis Plan ---")
    print("1. The key parent-daughter pair is Ba-140 -> La-140.")
    print("2. Assume Cherenkov counter only measures the daughter, La-140.")
    print("3. Use the two measurements to find the state of the system (A_Ba and A_La) at the first measurement.")
    print("4. Calculate the time (t1) required to reach this state from a pure Ba-140 sample.")
    print("\n--- Constants ---")
    print(f"Half-life of Ba-140: {T_Ba} days")
    print(f"Half-life of La-140: {T_La_hours} hours ({T_La:.4f} days)")
    print(f"Decay constant of Ba-140 (λ_Ba): {lambda_Ba:.4f} days^-1")
    print(f"Decay constant of La-140 (λ_La): {lambda_La:.4f} days^-1")
    print(f"Time between measurements (Δt): {delta_t} days")
    
    print("\n--- Calculation Steps ---")
    print(f"Activity of La-140 at measurement 1: {A_La_1} kBq/mL")
    print(f"Activity of La-140 at measurement 2: {A_La_2} kBq/mL")
    print("\nStep A: Find the activity of Ba-140 at the time of the first measurement (A_Ba_1).")
    print(f"Equation: A_La_2 = A_La_1 * e^(-λ_La*Δt) + A_Ba_1 * [λ_La/(λ_La-λ_Ba)] * [e^(-λ_Ba*Δt) - e^(-λ_La*Δt)]")
    print(f"Plugging in values: {A_La_2} = {A_La_1} * {exp_decay_La:.4f} + A_Ba_1 * {K:.4f} * ({exp_decay_Ba:.4f} - {exp_decay_La:.4f})")
    print(f"Solving for A_Ba_1 gives: {A_Ba_1:.4f} kBq/mL")

    print("\nStep B: Find the time (t1) from separation to the first measurement.")
    print(f"Using the activity ratio: A_La_1 / A_Ba_1 = {activity_ratio:.4f}")
    print("Equation: A_La_1 / A_Ba_1 = K * (1 - e^(-(λ_La - λ_Ba) * t1))")
    print(f"Plugging in values: {activity_ratio:.4f} = {K:.4f} * (1 - e^(-({lambda_La:.4f} - {lambda_Ba:.4f}) * t1))")
    print(f"Solving for t1 gives: {t1:.4f} days")

    # Final answer
    print("\n--- Conclusion ---")
    print(f"The approximate time between the sample separation and the first analysis is {t1:.2f} days.")

solve_radiochemistry_problem()
>>> 1.03