import math

def num_special_points(N):
    """
    Calculates the maximum possible number of special points for N planes
    in R^10, based on the argument that c=9.
    
    The formula is the binomial coefficient "N choose 9", C(N, 9).
    """
    k = 9
    if N < k:
        print(f"With N={N} planes, we cannot form an intersection of {k} planes.")
        print("Number of special points is 0.")
        return 0

    # Using math.comb for direct computation
    result = math.comb(N, k)
    
    print(f"For N = {N} planes, the maximum number of special points can be of the order O(N^c) where c = {k}.")
    print(f"This maximum is achieved in configurations where any {k} planes' direction vectors form a minimal spanning set.")
    print(f"The number of such points is bounded by C({N}, {k}).")
    print(f"C(N, k) = N! / (k! * (N-k)!)")
    
    # Let's show the numbers in the equation for a smaller example
    if N == 100:
      N_val = 100
      k_val = 9
      n_fact = math.factorial(N_val)
      k_fact = math.factorial(k_val)
      nmk_fact = math.factorial(N_val - k_val)
      print(f"For N={N_val}, k={k_val}:")
      # Using strings to represent potentially large numbers of factorials
      # To avoid showing gigantic numbers, we just format the expression.
      print(f"C({N_val}, {k_val}) = {N_val}! / ({k_val}! * ({N_val-k_val})!) = {N_val}! / ({k_val}! * {N_val-k_val}!)")
    
    print(f"Result: {result}")
    return result

# Example usage with N=100
num_special_points(100)