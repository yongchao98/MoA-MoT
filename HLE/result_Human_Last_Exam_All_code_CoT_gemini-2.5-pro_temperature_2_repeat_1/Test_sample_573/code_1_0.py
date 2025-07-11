import math

def calculate_simplices():
    """
    Calculates the number of n-simplices for N=200, k=13, and n from 0 to 5.
    
    The number of n-simplices in N_•(Z_N)_{k/} is the number of sequences
    (i_0, i_1, ..., i_n) such that k <= i_0 <= i_1 <= ... <= i_n <= N.
    This is a stars-and-bars problem, and the formula is C(N - k + n + 1, n + 1).
    """
    N = 200
    k = 13
    
    print(f"Calculating the number of n-simplices for N = {N}, k = {k}, and n <= 5.")
    
    results = []
    for n in range(6):
        # Upper argument of the binomial coefficient
        a = N - k + n + 1
        # Lower argument of the binomial coefficient
        b = n + 1
        
        result = math.comb(a, b)
        results.append(result)
        
        print(f"For n = {n}: Number of simplices = C({N} - {k} + {n} + 1, {n} + 1) = C({a}, {b}) = {result}")

    # The final answer in the requested format will be printed after the code execution.
    # The return value can be used to construct that answer.
    return results

if __name__ == "__main__":
    calculate_simplices()