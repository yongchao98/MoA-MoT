import math

def solve_ur_sequence():
    """
    Calculates and prints the sequence {u_r(3), u_r(4), ..., u_r(12)}.
    
    The value u_r(n) is the minimal order of the Picard-Fuchs differential
    equation for the real-analytic periods of the potential V(q) = (1/2)*(q^2 - q^n).
    The formula for n >= 3 is u_r(n) = floor((n + 1) / 2).
    """
    
    results = []
    print("The values of u_r(n) for n = 3 to 12 are calculated as follows:")
    
    # Loop through the required range of n
    for n in range(3, 13):
        # Apply the formula
        order = math.floor((n + 1) / 2)
        results.append(order)
        # Print each calculation as an "equation"
        print(f"u_r({n}) = floor(({n}+1)/2) = {order}")
        
    print("\nThe final sequence {u_r(3), u_r(4), ..., u_r(12)} is:")
    print(results)

solve_ur_sequence()