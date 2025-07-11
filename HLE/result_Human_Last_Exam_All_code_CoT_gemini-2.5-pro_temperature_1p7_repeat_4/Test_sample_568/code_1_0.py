import math

# --- Parameters ---
# Epsilon: the maximum allowed difference in concentration.
epsilon = 0.01
# Delta: the fraction of liquid exchanged in each step.
delta = 0.1

# --- Calculation ---
# This script calculates the number of steps 'n' required for the concentration
# of red liquid in the two cups to differ by at most epsilon.

# The derived formula is n = ceil(log(epsilon) / log(|1 - 2*delta|)).
# We must handle the edge case where delta = 0.5 separately to avoid a math error,
# as log(0) is undefined.
if not (0 < epsilon < 1 and 0 < delta < 1):
    print("Error: Epsilon and delta must be strictly between 0 and 1.")
    result = None
elif delta == 0.5:
    # If delta is 0.5, perfect mixing occurs in a single step.
    result = 1
else:
    # Applying the derived formula
    val_epsilon = epsilon
    val_delta_term = math.fabs(1 - 2 * delta)

    log_epsilon = math.log(val_epsilon)
    log_delta_term = math.log(val_delta_term)

    # The number of steps (as a float)
    n_float = log_epsilon / log_delta_term

    # Since the number of steps must be an integer, we take the ceiling.
    result = math.ceil(n_float)

# --- Output ---
print("The number of steps 'n' is determined by the formula:")
print("  n = ceil(ln(epsilon) / ln(|1 - 2*delta|))")
print("\nGiven the parameters:")
print(f"  epsilon = {epsilon}")
print(f"  delta   = {delta}")
print("\nThe final equation with substituted values is:")

if result is not None:
    if delta == 0.5:
        print("  n = 1 (This is a special case as delta=0.5 leads to perfect mixing in one step).")
    else:
        # Fulfilling the requirement to output each number in the final equation.
        print(f"  n = ceil(ln({val_epsilon}) / ln(|1 - 2*{delta}|))")
        print(f"  n = ceil({log_epsilon:.4f} / ln({val_delta_term:.4f}))")
        print(f"  n = ceil({log_epsilon:.4f} / {log_delta_term:.4f})")
        print(f"  n = ceil({n_float:.4f})")
        print(f"  n = {result}")
else:
    print("  Could not compute result due to invalid parameters.")
