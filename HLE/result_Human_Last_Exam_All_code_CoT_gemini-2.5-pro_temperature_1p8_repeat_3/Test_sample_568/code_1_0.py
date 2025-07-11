import math

def calculate_mixing_time(epsilon, delta):
    """
    Calculates the number of repetitions (n) needed for the cups to be "the same".

    Args:
        epsilon (float): The maximum allowed difference in liquid composition, must be in (0, 1).
        delta (float): The fraction of liquid moved in each step, must be in (0, 1).
    """
    # Validate inputs
    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon and delta must be strictly between 0 and 1.")
        return

    # Explain the model and formula
    print("Step 1: The recurrence relation for the fraction of red liquid in cup A (r_a) is:")
    print("r_a(n+1) = (1 - 2*\u03B4) * r_a(n) + \u03B4")
    print("\nStep 2: The solution to this is:")
    print("r_a(n) = 0.5 + 0.5 * (1 - 2*\u03B4)\u207F")
    print("\nStep 3: The 'sameness' condition |r_a(n) - r_b(n)| \u2264 \u03B5 simplifies to:")
    print("|1 - 2*\u03B4|\u207F \u2264 \u03B5")
    print("\nStep 4: Solving for n gives:")
    print("n \u2265 log(\u03B5) / log(|1 - 2*\u03B4|)")
    print("So, we need to calculate n = ceil(log(\u03B5) / log(|1 - 2*\u03B4|))\n")


    print("-------------------- CALCULATION --------------------")
    print(f"Given parameters: \u03B5 = {epsilon}, \u03B4 = {delta}")
    
    # Special case: If delta is 0.5, perfect mixing occurs in one step.
    if delta == 0.5:
        n = 1
        print("\nSpecial case: Since \u03B4 = 0.5, the liquids mix perfectly in a single step.")
        final_n = n
    else:
        # General case calculation
        abs_val = abs(1 - 2 * delta)
        log_epsilon = math.log(epsilon)
        log_abs_val = math.log(abs_val)

        # The result of the division
        ratio = log_epsilon / log_abs_val
        
        # n is the ceiling of this ratio
        final_n = math.ceil(ratio)

        print("\nSubstituting the values into the formula:")
        print(f"n \u2265 log({epsilon}) / log(|1 - 2*{delta}|)")
        print(f"n \u2265 {log_epsilon:.4f} / log({abs_val:.4f})")
        print(f"n \u2265 {log_epsilon:.4f} / {log_abs_val:.4f}")
        print(f"n \u2265 {ratio:.4f}")
        print(f"\nThe smallest integer n that satisfies this is ceil({ratio:.4f}).")

    print("\n-------------------- FINAL ANSWER --------------------")
    print(f"The procedure must be repeated {final_n} times.")

if __name__ == '__main__':
    # --- Parameters ---
    # You can change these values to match your specific problem.
    epsilon_param = 0.01  # example value for ε
    delta_param = 0.1     # example value for δ
    # --- End of Parameters ---
    
    calculate_mixing_time(epsilon_param, delta_param)