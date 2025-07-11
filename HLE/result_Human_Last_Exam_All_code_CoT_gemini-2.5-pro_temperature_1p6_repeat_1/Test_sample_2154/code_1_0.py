import math

def solve_picard_fuchs_order():
    """
    Calculates the minimal order of the Picard-Fuchs differential equation u_r(n)
    for the Hamiltonian V(q) = (q^2 - q^n)/2 for n from 3 to 12.

    The order is given by the formula u_r(n) = floor((n-1)/2).
    """
    print("The values for u_r(n) are calculated using the formula: u_r(n) = floor((n-1)/2)")
    print("-" * 70)
    
    results = []
    for n in range(3, 13):
        # Calculate the order using the formula
        numerator = n - 1
        denominator = 2
        # In Python, integer division // is equivalent to floor division for positive results.
        u_r_n = numerator // denominator
        
        # Store the result
        results.append(u_r_n)
        
        # As requested, output each number in the final equation for each step.
        print(f"For n = {n}: u_r({n}) = floor(({n}-1)/2) = floor({numerator}/2) = {u_r_n}")

    # The problem asks to find the set {u_r(3), u_r(4), ..., u_r(12)}.
    # The script has printed each individual value, now let's print the final set.
    print("-" * 70)
    print(f"The complete set is: {results}")

solve_picard_fuchs_order()