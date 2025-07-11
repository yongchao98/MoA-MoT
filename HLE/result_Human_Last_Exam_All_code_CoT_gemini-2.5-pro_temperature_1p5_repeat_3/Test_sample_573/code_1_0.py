import math

def solve_simplices_count():
    """
    Calculates the number of n-simplices for the given mathematical object.
    
    An n-simplex in N_.(Z_N)_k/ corresponds to a sequence of integers
    (c_0, c_1, ..., c_n) such that -N <= c_0 <= c_1 <= ... <= c_n <= k.
    
    The number of such sequences is given by the combinatorial formula:
    C( (k - (-N) + 1) + (n+1) - 1, n+1 ) = C(k + N + n + 1, n + 1)
    """
    N = 200
    k = 13

    print(f"Calculating the number of n-simplices for N={N} and k={k}, for n from 0 to 5.")
    print("-" * 30)
    
    results = []
    for n in range(6):
        # Upper number in the combination C(total, choose)
        total = k + N + n + 1
        # Lower number in the combination C(total, choose)
        choose = n + 1
        
        # Calculate the number of simplices using math.comb
        num_simplices = math.comb(total, choose)
        results.append(num_simplices)
        
        # Print the result for each n, showing the formula
        print(f"For n = {n}, the number of simplices is C({total}, {choose}) = {num_simplices}")
        
    # Final answer in the required format
    answer_str = ",".join(map(str, results))
    print(f"\n<<<__{answer_str}__>>>")

solve_simplices_count()