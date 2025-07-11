import math

def calculate_simplices(N, k, n_max):
    """
    Calculates the number of n-simplices for n from 0 to n_max.

    An n-simplex in N.(Z_N)_{k/} corresponds to a sequence of integers
    (j_0, j_1, ..., j_n) such that k <= j_0 <= j_1 <= ... <= j_n <= N.

    This is a combination with repetition problem. We are choosing n+1 numbers
    from the set {k, k+1, ..., N}, which has (N - k + 1) elements.
    The number of such sequences is given by the binomial coefficient:
    C((N - k + 1) + (n + 1) - 1, n + 1) = C(N - k + n + 1, n + 1)
    """
    
    print(f"Calculating the number of n-simplices for N={N}, k={k}, and n up to {n_max}.\n")
    print("The number of n-simplices is given by the formula: C(N - k + n + 1, n + 1)")
    print("-" * 60)

    # Number of choices for each integer in the sequence
    num_choices = N - k + 1
    
    results = []
    for n in range(n_max + 1):
        # Number of integers in the sequence
        num_to_choose = n + 1
        
        # Binomial coefficient C(p + m - 1, m)
        comb_n = num_choices + num_to_choose - 1
        comb_k = num_to_choose
        
        # In terms of N, k, n:
        # comb_n = (N - k + 1) + (n + 1) - 1 = N - k + n + 1
        # comb_k = n + 1
        
        result = math.comb(comb_n, comb_k)
        results.append(result)
        
        print(f"For n = {n}:")
        print(f"Number of simplices = C({N} - {k} + {n} + 1, {n} + 1) = C({comb_n}, {comb_k}) = {result}")

    return results

if __name__ == '__main__':
    N = 200
    k = 13
    n_max = 5
    
    calculate_simplices(N, k, n_max)
