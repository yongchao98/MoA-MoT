import math

def calculate_simplices(N, k, n_max):
    """
    Calculates and prints the number of n-simplices for n up to n_max.

    The number of n-simplices in the nerve of the overcategory Z_N over k
    is the number of non-decreasing sequences of length n+1 of integers
    in the range [k, N]. This is a combination with repetition problem.

    The formula is C(m + r - 1, r), where m is the number of items to choose from
    and r is the number of items being chosen.
    Here, m = N - k + 1 (size of the set {k, ..., N})
    and r = n + 1 (length of the sequence x_0, ..., x_n).
    So the formula is C((N - k + 1) + (n + 1) - 1, n + 1) = C(N - k + n + 1, n + 1).
    """

    print(f"Calculating the number of n-simplices for N = {N}, k = {k}, and n <= {n_max}.\n")
    
    for n in range(n_max + 1):
        # Parameters for the binomial coefficient C(a, b)
        a = N - k + n + 1
        b = n + 1
        
        # Calculate the result
        try:
            result = math.comb(a, b)
            print(f"For n = {n}:")
            # The problem asks to output each number in the final equation.
            # Here is the full equation with all numbers shown.
            equation = f"C({N} - {k} + {n} + 1, {n} + 1)"
            simplified_equation = f"C({a}, {b})"
            print(f"  Number of {n}-simplices = {equation} = {simplified_equation} = {result}")
        except ValueError:
            print(f"Cannot calculate C({a}, {b}) for n = {n}. Invalid input.")

if __name__ == '__main__':
    N_val = 200
    k_val = 13
    n_max_val = 5
    calculate_simplices(N_val, k_val, n_max_val)
