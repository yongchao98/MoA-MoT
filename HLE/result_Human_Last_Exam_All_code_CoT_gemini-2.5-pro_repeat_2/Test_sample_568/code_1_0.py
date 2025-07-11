import math

# Parameters for the problem. You can change these values.
epsilon = 0.01
delta = 0.2

# The problem states epsilon, delta are in (0, 1). We assume valid inputs.
# The derived formula is t = ceil(log(epsilon / 2) / log(|1 - 2*delta|))

# Handle the special case of delta = 0.5, which results in perfect mixing in 1 step.
if delta == 0.5:
    t_final = 1
    # Print the equation steps for this simple case
    print("For delta = 0.5, the difference becomes 0 after 1 step.")
    print("t = 1")
else:
    # --- Step-by-step evaluation of the formula ---
    
    # Values for the equation
    val_eps_div_2 = epsilon / 2
    val_1_minus_2_delta_abs = abs(1 - 2 * delta)
    
    # Calculate numerator and denominator
    numerator = math.log(val_eps_div_2)
    denominator = math.log(val_1_minus_2_delta_abs)
    
    # Calculate the floating point result before ceiling
    t_float = numerator / denominator
    
    # Calculate the final integer result
    t_final = math.ceil(t_float)

    # --- Outputting each step of the calculation as requested ---
    
    # Print the equation with initial parameters
    print(f"t = ceil(log({epsilon} / 2) / log(|1 - 2 * {delta}|))")
    
    # Print the equation after initial arithmetic
    print(f"t = ceil(log({val_eps_div_2}) / log({val_1_minus_2_delta_abs}))")

    # Print the equation with logarithm values
    print(f"t = ceil({numerator} / {denominator})")

    # Print the equation after division
    print(f"t = ceil({t_float})")

    # Print the final result
    print(f"t = {t_final}")
