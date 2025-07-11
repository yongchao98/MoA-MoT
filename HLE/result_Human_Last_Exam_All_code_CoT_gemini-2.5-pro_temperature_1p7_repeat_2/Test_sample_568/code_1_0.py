import math

# --- Parameters ---
# You can change these values to see the result for different parameters.
epsilon = 0.01
delta = 0.1

# --- Calculation ---

print(f"Calculating the number of repetitions for epsilon = {epsilon} and delta = {delta}\n")

# Validate parameters
if not (0 < epsilon < 1):
    print("Error: Epsilon (epsilon) must be in the range (0, 1).")
elif not (0 < delta < 1):
    print("Error: Delta (delta) must be in the range (0, 1).")

# Handle the special case where delta = 0.5
elif delta == 0.5:
    n = 1
    print("For delta = 0.5, the liquids mix perfectly in one step.")
    print(f"The number of repetitions n = {n}")

# General case for delta in (0, 1) but not 0.5
else:
    print("The formula for the number of steps n is:")
    print("n = ceil(ln(epsilon) / ln(|1 - 2*delta|))\n")
    print("--- Calculation Breakdown ---")
    
    # Calculate the numerator
    log_epsilon = math.log(epsilon)
    print(f"ln(epsilon) = ln({epsilon}) = {log_epsilon}")
    
    # Calculate the denominator
    abs_val = abs(1 - 2 * delta)
    print(f"|1 - 2*delta| = |1 - 2*{delta}| = |{1 - 2 * delta}| = {abs_val}")
    
    log_abs_val = math.log(abs_val)
    print(f"ln(|1 - 2*delta|) = ln({abs_val}) = {log_abs_val}")
    
    # Final calculation for n
    result_float = log_epsilon / log_abs_val
    n = math.ceil(result_float)
    
    print("\n--- Final Result ---")
    print(f"n = ceil({log_epsilon} / {log_abs_val})")
    print(f"n = ceil({result_float})")
    print(f"The number of repetitions n = {n}")
