import math

def calculate_simplices():
    """
    Calculates the number of n-simplices for n from 0 to 5,
    with N=200 and k=13, based on the derived combinatorial formula.
    """
    N = 200
    k = 13

    # The number of integers in the range [k, N]
    num_choices = N - k + 1
    
    print(f"For N = {N} and k = {k}:")
    print(f"The number of available integers is N - k + 1 = {N} - {k} + 1 = {num_choices}.")
    print("The number of n-simplices is calculated by the formula C(num_choices + n, n + 1).")
    print("-" * 50)
    
    results = []
    for n in range(6):
        # We need to choose n+1 elements for an n-simplex
        # from 'num_choices' available integers with repetition.
        # This is C(num_choices + (n+1) - 1, n+1) = C(num_choices + n, n+1)
        
        n_comb = num_choices + n
        k_comb = n + 1
        
        result = math.comb(n_comb, k_comb)
        results.append(result)
        
        print(f"For n = {n}:")
        print(f"  Number of simplices = C({num_choices} + {n}, {n} + 1) = C({n_comb}, {k_comb}) = {result}")

calculate_simplices()