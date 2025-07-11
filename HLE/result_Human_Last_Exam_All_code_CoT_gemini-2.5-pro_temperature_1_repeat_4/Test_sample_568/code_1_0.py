import math

# --- Parameters ---
# Epsilon (epsilon): The maximum allowed difference in concentration.
# Must be in the range (0, 1).
epsilon = 0.01

# Delta (delta): The fraction of liquid moved in each step.
# Must be in the range (0, 1).
delta = 0.25
# ------------------

print(f"Goal: Find the number of steps 'n' for epsilon = {epsilon} and delta = {delta}.")
print("-" * 20)

# Handle the special case where delta is 0.5
if delta == 0.5:
    # If delta is 0.5, perfect mixing occurs in one step.
    # The difference becomes 0, which is always <= epsilon.
    n_repetitions = 1
    print("For delta = 0.5, perfect mixing is achieved in a single step.")
    print(f"The number of repetitions n is: {n_repetitions}")

else:
    # Check for invalid delta values that would make the calculation impossible
    if delta <= 0 or delta >= 1:
        print("Error: delta must be between 0 and 1.")
    else:
        # The number of repetitions 'n' must satisfy the inequality:
        # |1 - 2*delta|^n <= epsilon
        # which is equivalent to:
        # n >= log(epsilon) / log(|1 - 2*delta|)
        
        # Calculate the terms for the formula
        log_epsilon = math.log(epsilon)
        abs_val = abs(1 - 2 * delta)
        log_abs_val = math.log(abs_val)
        
        # Calculate the lower bound for n
        n_float = log_epsilon / log_abs_val
        
        # Since n must be an integer, we take the ceiling of the result
        n_repetitions = math.ceil(n_float)

        print(f"The number of repetitions 'n' must satisfy the inequality: |1 - 2 * {delta}|^n <= {epsilon}")
        print(f"This is equivalent to: n >= log({epsilon}) / log(|1 - 2 * {delta}|)")
        print(f"n >= log({epsilon}) / log({abs_val})")
        print(f"n >= {log_epsilon} / {log_abs_val}")
        print(f"n >= {n_float}")
        print("\nSince n must be an integer, we take the ceiling.")
        print(f"The smallest integer n is: {n_repetitions}")

print("-" * 20)
print(f"Final answer: The procedure needs to be repeated {n_repetitions} times.")