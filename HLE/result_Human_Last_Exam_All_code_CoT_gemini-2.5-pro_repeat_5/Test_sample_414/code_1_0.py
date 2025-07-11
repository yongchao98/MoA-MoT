import math

def solve_radiochemistry_problem():
    """
    This script solves for the time between sample irradiation and the first analysis
    based on the principles of radioactive parent-daughter decay.
    """
    
    # --- Step 1: Define constants based on the problem description ---
    
    # Half-life of the parent nuclide, Ba-140, in days
    T_half_parent = 12.75
    # Half-life of the daughter nuclide, La-140, in days
    T_half_daughter = 1.678
    
    # Measured activities in kBq/mL
    A1 = 1.4  # Activity at the unknown time 't'
    A2 = 2.1  # Activity at time 't + 14 days'
    
    # Time interval between the two measurements in days
    delta_t = 14.0

    # --- Step 2: Calculate decay constants ---
    # The decay constant lambda is calculated as ln(2) / T_half
    lambda_parent = math.log(2) / T_half_parent
    lambda_daughter = math.log(2) / T_half_daughter

    # --- Step 3: Solve for the unknown time 't' ---
    # We set up a ratio of the activities at t and t + delta_t.
    # A2 / A1 = (exp(-lp*t)*exp(-lp*dt) - exp(-ld*t)*exp(-ld*dt)) / (exp(-lp*t) - exp(-ld*t))
    # After rearranging, we can solve for t:
    # t = ln( (R - exp(-ld*dt)) / (R - exp(-lp*dt)) ) / (ld - lp)
    # where R is the activity ratio A2/A1.

    activity_ratio = A2 / A1
    
    exp_term_parent = math.exp(-lambda_parent * delta_t)
    exp_term_daughter = math.exp(-lambda_daughter * delta_t)

    # Numerator of the argument for the natural logarithm
    log_arg_numerator = activity_ratio - exp_term_daughter
    
    # Denominator of the argument for the natural logarithm
    log_arg_denominator = activity_ratio - exp_term_parent
    
    log_argument = log_arg_numerator / log_arg_denominator
    
    # Difference in decay constants
    lambda_diff = lambda_daughter - lambda_parent

    # Final calculation for time 't'
    time_t = math.log(log_argument) / lambda_diff

    # --- Step 4: Print the results clearly ---
    print("Based on the Ba-140/La-140 decay chain, we solve for the time 't'.")
    print("The governing equation is derived from the ratio of activities:")
    print("t = ln((A2/A1 - exp(-lambda_La * delta_t)) / (A2/A1 - exp(-lambda_Ba * delta_t))) / (lambda_La - lambda_Ba)")
    print("\nSubstituting the calculated numerical values into the equation:")
    
    # This print statement shows the final equation with all numbers plugged in.
    final_equation_str = (
        f"t = ln(({activity_ratio:.3f} - {exp_term_daughter:.5f}) / "
        f"({activity_ratio:.3f} - {exp_term_parent:.5f})) / "
        f"({lambda_daughter:.5f} - {lambda_parent:.5f})"
    )
    print(final_equation_str)
    
    print("\nWhich simplifies to:")
    final_calc_str = f"t = ln({log_argument:.4f}) / {lambda_diff:.5f}"
    print(final_calc_str)
    
    print(f"t = {math.log(log_argument):.4f} / {lambda_diff:.5f}")

    print("\n--- Final Answer ---")
    print(f"The approximate time between irradiation and the first analysis is: {time_t:.2f} days")

solve_radiochemistry_problem()