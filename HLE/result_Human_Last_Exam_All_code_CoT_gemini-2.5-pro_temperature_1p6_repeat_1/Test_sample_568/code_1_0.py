import math

# You can change these parameters to test different scenarios.
epsilon = 0.01
delta = 0.1

print(f"Calculating the number of steps 'n' for epsilon = {epsilon} and delta = {delta}\n")

# Handle the edge case where delta = 0.5
if delta == 0.5:
    # If delta is 0.5, the cups mix perfectly in one step. The concentration
    # difference becomes 0. The general formula is undefined in this case
    # because it would involve log(0).
    n = 1
    print("For delta = 0.5, the mixing is completed in a single step.")
    print(f"The required number of steps is {n}.")

else:
    # This is the general case for delta in (0, 1) but not equal to 0.5.
    # The formula to find the number of steps n is:
    # n = ceil(log(epsilon) / log(|1 - 2*delta|))

    # Calculate the components of the formula
    val_1_minus_2_delta = 1 - 2 * delta
    base = abs(val_1_minus_2_delta)
    
    log_epsilon = math.log(epsilon)
    log_base = math.log(base)
    
    n_float = log_epsilon / log_base
    n = math.ceil(n_float)

    # Print the final equation and the intermediate values, as requested.
    # The requirement is to "output each number in the final equation".
    print("The number of steps 'n' is calculated using the formula:")
    print("n = ceil(log(epsilon) / log(|1 - 2*delta|))\n")
    
    print("Substituting the given values:")
    print(f"n = ceil(log({epsilon}) / log(|1 - 2*{delta}|))")
    print(f"n = ceil(log({epsilon}) / log({base:.4f}))")
    print(f"n = ceil({log_epsilon:.4f} / {log_base:.4f})")
    print(f"n = ceil({n_float:.4f})")
    print(f"n = {n}\n")
    
    print(f"The final answer is that the procedure needs to be repeated {n} times.")
