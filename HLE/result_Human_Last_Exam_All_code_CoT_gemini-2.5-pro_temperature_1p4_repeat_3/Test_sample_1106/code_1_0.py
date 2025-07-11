import math

def calculate_min_max_points(n, d):
    """
    Calculates the minimized maximum number of points in any closed hyper-hemisphere.
    
    Args:
        n (int): The number of points.
        d (int): The dimension of the hypersphere.
    """
    
    # Based on known results in combinatorial geometry, the formula depends
    # on the parity of n and d.
    
    # Case 1: n is even and d is odd
    if n % 2 == 0 and d % 2 != 0:
        val1 = n // 2
        val2 = d // 2
        result = val1 + val2
        print(f"For N={n} (even) and d={d} (odd), the formula is N/2 + floor(d/2).")
        print(f"Calculation: {n}/2 + floor({d}/2) = {val1} + {val2} = {result}")
    # Case 2: All other combinations of parity for n and d
    else:
        val1 = n // 2
        val2 = math.ceil(d / 2)
        result = val1 + val2
        
        parity_n = "odd" if n % 2 != 0 else "even"
        parity_d = "odd" if d % 2 != 0 else "even"

        print(f"For N={n} ({parity_n}) and d={d} ({parity_d}), the applicable formula is floor(N/2) + ceil(d/2).")
        print(f"Calculation: floor({n}/2) + ceil({d}/2) = {val1} + {val2} = {result}")

    print("\nThe largest number of points that can be achieved in the minimized-maximum arrangement is:")
    print(result)
    return result

# Given values from the problem
N = 15
D = 8

calculate_min_max_points(N, D)
