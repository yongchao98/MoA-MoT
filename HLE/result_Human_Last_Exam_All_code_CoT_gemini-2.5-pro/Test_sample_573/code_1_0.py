import math

def calculate_simplices(n, N, k):
    """
    Calculates the number of n-simplices in N_.(Z_N)_{k/}.
    The formula is C(N - k + n + 1, n + 1).
    """
    # The number of objects to choose from is N - k + 1
    # The number of items to choose is n + 1
    a = N - k + n + 1
    b = n + 1
    
    # Using math.comb for combinations C(a, b)
    result = math.comb(a, b)
    
    return a, b, result

def main():
    """
    Main function to solve the problem for the given parameters.
    """
    N = 200
    k = 13
    n_values = range(6) # for n = 0, 1, 2, 3, 4, 5

    print(f"Calculating the number of n-simplices for N={N}, k={k}, and n <= 5.\n")
    print("The number of n-simplices is given by the formula C(N - k + n + 1, n + 1).")
    print("-" * 70)
    
    results = []
    for n in n_values:
        a, b, result = calculate_simplices(n, N, k)
        print(f"For n = {n}, the number of simplices is C({N} - {k} + {n} + 1, {n} + 1) = C({a}, {b}) = {result}")
        results.append(result)

if __name__ == "__main__":
    main()
