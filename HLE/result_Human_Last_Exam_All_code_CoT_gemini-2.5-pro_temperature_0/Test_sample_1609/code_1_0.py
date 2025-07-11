import math

def solve_max_prime_implicants():
    """
    Calculates a(4), the maximal number of prime implicants of a Boolean function of 4 variables.
    """
    n = 4

    print(f"To find a(n), the maximal number of prime implicants for a Boolean function of n variables, we use a known result from combinatorics.")
    print("The maximum value is given by the largest binomial coefficient C(n, k) for a fixed n.")
    print("This coefficient is maximized when k is as close to n/2 as possible, which is k = floor(n/2).")
    print(f"So, the formula is a(n) = C(n, floor(n/2)).")

    # Step 1: Set n
    print(f"\nFor this problem, n = {n}.")

    # Step 2: Calculate k
    k = math.floor(n / 2)
    print(f"The value of k is floor({n}/2) = {k}.")

    # Step 3: Calculate the result using the formula a(n) = C(n, k)
    result = math.comb(n, k)

    # Step 4: Print the final equation with all numbers
    print("\nThe final calculation is:")
    print(f"a({n}) = C({n}, {k}) = {result}")

solve_max_prime_implicants()