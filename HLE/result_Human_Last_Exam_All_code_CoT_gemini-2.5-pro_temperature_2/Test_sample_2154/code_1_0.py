import math

def solve_u_r_sequence():
    """
    Calculates and prints the values of u_r(n) for n from 3 to 12.
    u_r(n) is the minimal order of the Picard-Fuchs differential equation for
    the periods of the Hamiltonian V(q) = 1/2 * (q^2 - q^n).
    The formula, based on the genus of the associated hyperelliptic curve, is
    u_r(n) = 2 * floor((n-1)/2).
    """
    
    results = []
    print("Calculating the values of u_r(n) for n = 3 to 12:")
    
    for n in range(3, 13):
        # Perform calculations step-by-step for clarity
        numerator = n - 1
        division_result = numerator / 2
        floor_value = math.floor(division_result)
        final_value = 2 * floor_value
        results.append(final_value)
        
        # Output each number in the calculation as requested
        print(f"u_r({n}) = 2 * floor(({n}-1)/2) = 2 * floor({numerator}/2) = 2 * floor({division_result}) = 2 * {floor_value} = {final_value}")

    print(f"\nThe final list of values for {{u_r(3), u_r(4), ..., u_r(12)}} is:\n{results}")

solve_u_r_sequence()