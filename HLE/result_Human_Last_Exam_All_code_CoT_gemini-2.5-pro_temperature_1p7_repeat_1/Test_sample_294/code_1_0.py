import math

def solve_k_matching_complexity():
    """
    This function calculates the maximum integer k for which counting k-matchings
    is known to be solvable in subcubic time, based on established algorithms.

    The time complexity of the best known algorithm is O(n^(omega * ceil(k/3))),
    where n is the number of vertices and omega is the matrix multiplication exponent.
    We need to find the largest k for which this exponent is less than 3.
    """
    
    # The current best upper bound for the matrix multiplication exponent.
    omega = 2.3728596 
    
    # The cubic time barrier we want to stay under.
    cubic_exponent = 3.0
    
    max_k = 0
    k = 1
    
    print("Analyzing the complexity of counting k-matchings: O(n^(omega * ceil(k/3)))")
    print(f"Using omega = {omega}")
    print(f"Seeking the largest integer k such that exponent < {cubic_exponent}")
    print("-" * 60)
    
    while True:
        # Calculate ceil(k/3)
        ceil_k_div_3 = math.ceil(k / 3)
        
        # Calculate the overall exponent
        exponent = omega * ceil_k_div_3
        
        is_subcubic = exponent < cubic_exponent
        
        print(f"For k = {k}:")
        # I was instructed to output each number in the final equation.
        # Equation: omega * ceil(k/3) = exponent
        print(f"  Calculation: {omega:.5f} * ceil({k}/3) = {omega:.5f} * {ceil_k_div_3} = {exponent:.5f}")

        if is_subcubic:
            print(f"  Result: Exponent {exponent:.5f} is < {cubic_exponent}. The algorithm is subcubic.")
            max_k = k
            k += 1
        else:
            print(f"  Result: Exponent {exponent:.5f} is >= {cubic_exponent}. The algorithm is not subcubic.")
            print("-" * 60)
            print(f"The analysis stops. The maximum k for which the known algorithm is subcubic is {max_k}.")
            break
            
solve_k_matching_complexity()
