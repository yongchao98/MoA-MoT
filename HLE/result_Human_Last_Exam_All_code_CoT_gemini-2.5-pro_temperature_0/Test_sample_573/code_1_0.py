import math

def count_simplices():
    """
    Calculates the number of n-simplices for the nerve of the overcategory Z_N over k.
    
    The number of n-simplices in N_•(Z_N)_k/ is the number of sequences
    k <= k_0 <= k_1 <= ... <= k_n <= N.
    This is a combinations with repetition problem, and the formula is:
    C( (N - k + 1) + (n + 1) - 1, n + 1) = C(N - k + n + 1, n + 1).
    """
    N = 200
    k = 13

    print(f"Calculating the number of n-simplices for N = {N} and k = {k}.")
    print("The formula is: C(N - k + n + 1, n + 1)")
    print("-" * 50)

    for n in range(6):
        # Parameters for the combination formula C(a, b)
        a = N - k + n + 1
        b = n + 1
        
        # Calculate the number of simplices using math.comb
        result = math.comb(a, b)
        
        # Print the detailed calculation for each n
        print(f"For n = {n}:")
        equation = f"Number of {n}-simplices = C({N} - {k} + {n} + 1, {n} + 1)"
        calculation = f"C({a}, {b})"
        print(f"{equation} = {calculation} = {result}")
        print()

count_simplices()