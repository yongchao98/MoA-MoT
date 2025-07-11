import math

# --- Parameters ---
# You can change these values to solve for different epsilon and delta.
epsilon = 0.01  # Epsilon: The maximum allowed difference (e.g., 0.01 for 1%)
delta = 0.1   # Delta: The fraction of liquid moved in each step (e.g., 0.1 for 10%)

# --- Plan ---
# The number of steps 'n' is calculated based on the formula derived from the mixing process.
# The general formula is: n = ceil(log(epsilon) / log(abs(1 - 2*delta)))
# A special case exists for delta = 0.5, where perfect mixing occurs in one step.

print(f"Goal: Find the number of steps 'n' for epsilon = {epsilon} and delta = {delta}\n")

# --- Calculation ---

# Handle the special case where delta = 0.5
if delta == 0.5:
    n = 1
    print("This is a special case where delta = 0.5.")
    print("The liquids become perfectly mixed in a single step.")
    print("\nFinal Answer:")
    print(f"n = {n}")

else:
    # Calculate using the formula for the general case.
    # Step 1: Calculate the numerator, log(epsilon)
    log_epsilon = math.log(epsilon)

    # Step 2: Calculate the base of the second logarithm, |1 - 2*delta|
    base_val = abs(1 - 2 * delta)

    # Step 3: Calculate the denominator, log(|1 - 2*delta|)
    log_base = math.log(base_val)

    # Step 4: Perform the division
    ratio = log_epsilon / log_base

    # Step 5: Take the ceiling to find the smallest integer n
    n = int(math.ceil(ratio))

    # --- Output the Calculation Steps ---
    print("The number of steps 'n' is the smallest integer satisfying: |1 - 2*delta|^n <= epsilon")
    print("This is solved by: n = ceil(log(epsilon) / log(abs(1 - 2*delta)))")
    print("-" * 30)
    
    # Show equation with substituted values
    print(f"n = ceil(log({epsilon}) / log(abs(1 - 2*{delta})))")
    
    # Show intermediate calculations
    print(f"n = ceil({log_epsilon:.4f} / log(abs({1 - 2*delta:.2f})))")
    print(f"n = ceil({log_epsilon:.4f} / log({base_val:.2f}))")
    print(f"n = ceil({log_epsilon:.4f} / {log_base:.4f})")
    print(f"n = ceil({ratio:.4f})")

    # Show final result
    print("\nFinal Answer:")
    print(f"n = {n}")
