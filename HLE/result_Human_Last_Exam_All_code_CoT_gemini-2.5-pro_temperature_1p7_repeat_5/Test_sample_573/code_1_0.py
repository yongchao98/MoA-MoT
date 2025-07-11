import math

def calculate_simplices():
    """
    Calculates the number of n-simplices for the given problem.

    The number of n-simplices in the nerve of the overcategory N_•(Z_N)_{k/}
    is the number of non-decreasing sequences (c_0, c_1, ..., c_n)
    such that k <= c_0 <= c_1 <= ... <= c_n <= N.

    This is a combination with repetition problem. The number of such sequences
    is given by C(m + r - 1, r), where m is the number of items to choose from
    and r is the number of choices.
    Here, m = (N - k + 1) is the number of objects in the overcategory,
    and r = n + 1 is the length of the simplex sequence.
    The formula is C((N - k + 1) + (n + 1) - 1, n + 1) = C(N - k + n + 1, n + 1).
    """
    N = 200
    k = 13
    n_max = 5

    # Number of objects in the overcategory (Z_N)_{k/}
    num_objects = N - k + 1

    print(f"Given parameters: N = {N}, k = {k}")
    print(f"The number of objects in the overcategory (Z_N)_{k/} is {N} - {k} + 1 = {num_objects}.")
    print("The number of n-simplices is calculated using the formula: C(num_objects + n, n + 1)\n")

    results = []
    for n in range(n_max + 1):
        # We are choosing n+1 numbers from num_objects with replacement.
        # This is C(num_objects + n, n + 1)
        m = num_objects + n
        r = n + 1
        
        num_simplices = math.comb(m, r)
        results.append(num_simplices)
        
        print(f"For n = {n}:")
        # Print the equation with all the numbers
        equation = f"C({num_objects} + {n}, {n} + 1) = C({m}, {r})"
        print(f"  The number of {n}-simplices is {equation} = {num_simplices:,}")

if __name__ == '__main__':
    calculate_simplices()
