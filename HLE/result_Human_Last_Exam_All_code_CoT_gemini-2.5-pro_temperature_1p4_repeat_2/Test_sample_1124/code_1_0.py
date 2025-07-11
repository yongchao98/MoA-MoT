import sys

def modified_logistic_map(x, r):
    """
    A modified logistic map with a stable fixed point.
    The equation is: X_n+1 = X_n - (X_n - 1) / R
    """
    # This check prevents division by zero, though R=0 is not relevant for this problem.
    if r == 0:
        # In this case, the equation is undefined. Return x to indicate no change.
        return x
    return x - (x - 1) / r

def find_equilibrium():
    """
    Demonstrates the equilibrium point of the modified map.
    """
    r = 3.57
    # Start with a value different from the equilibrium to show convergence
    x = 2.5 
    
    # Iterate a few times to converge to the fixed point
    for _ in range(20):
        x = modified_logistic_map(x, r)
        
    # The converged value is the equilibrium point
    equilibrium_point = x
    
    # The final equation requires outputting each number.
    # At equilibrium, X_n+1 = X_n = equilibrium_point.
    # We will show one step of the calculation at the equilibrium point.
    x_next = modified_logistic_map(equilibrium_point, r)
    
    print("The modified map is: X_n+1 = X_n - (X_n - 1) / R")
    print(f"\nFor R = {r}, the system gets a stable equilibrium point.")
    print("\nDemonstrating the equilibrium condition (X_n+1 = X_n):")
    # Using high precision for demonstration, though the result is exact.
    # The format is X_n+1 = calculation = result
    print(f"{x_next:.8f} = {equilibrium_point:.8f} - ({equilibrium_point:.8f} - 1) / {r}")
    
    # Returning the value for the final answer block
    return equilibrium_point

# Execute the function to print the results
final_answer = find_equilibrium()
# Suppress the final print of the object itself to stdout for cleaner output
# and prepare the value for the <<<>>> block.
# sys.stdout = open(os.devnull, 'w') # this is not needed if running in a standard env

# The value should be printed to the screen though
# The final answer required by the format.
# print(f"<<<{final_answer}>>>")