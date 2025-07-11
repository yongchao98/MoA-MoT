import math

def count_non_admissible(a, b):
    """
    Calculates the number of non-admissible integers k for given positive integers a and b.

    An integer k is "admissible" if there exist complex a x b matrices A_1,...,A_{ab} satisfying:
    1. Each A_i is nonzero.
    2. tr(A_i^dagger A_j) = 0 whenever i != j (orthogonality).
    3. Exactly k of the matrices A_i have rank 1.

    Args:
        a (int): The number of rows, a positive integer.
        b (int): The number of columns, a positive integer.
    """
    if not (isinstance(a, int) and a > 0 and isinstance(b, int) and b > 0):
        print(f"Error: a and b must be positive integers. Got a={a}, b={b}")
        return

    ab = a * b
    
    # Case 1: min(a, b) = 1
    # In this case, any nonzero a x b matrix must have rank 1.
    # Therefore, any orthogonal basis of nonzero matrices must consist of ab rank-1 matrices.
    # So, k=ab is the only admissible value.
    # The non-admissible values are 0, 1, ..., ab-1.
    if min(a, b) == 1:
        num_non_admissible = ab
        print(f"For a={a}, b={b}:")
        print(f"Since min(a, b) = 1, any nonzero matrix must have rank 1.")
        print(f"The only admissible value for k is {ab}.")
        print(f"The non-admissible integers are 0, 1, ..., {ab-1}.")
        print(f"Total number of non-admissible integers = {num_non_admissible}")
        print("-" * 20)

    # Case 2: min(a, b) >= 2
    # In this case, based on results from matrix theory, it's known that
    # the only non-admissible values for k are 1 and ab-1.
    else:
        num_non_admissible = 2
        non_admissible_k1 = 1
        non_admissible_k2 = ab - 1
        print(f"For a={a}, b={b}:")
        print(f"Since min(a, b) >= 2, the set of admissible values for k includes all integers")
        print(f"in the range [0, {ab}] except for two specific values.")
        print(f"The non-admissible integers are k = {non_admissible_k1} and k = {ab} - 1 = {non_admissible_k2}.")
        print(f"Total number of non-admissible integers = {num_non_admissible}")
        print("-" * 20)

def main():
    """
    Runs the calculation for a few example pairs of (a, b).
    """
    # Example 1: a=1, b=5 (Case 1)
    count_non_admissible(1, 5)

    # Example 2: a=3, b=1 (Case 1)
    count_non_admissible(3, 1)

    # Example 3: a=2, b=2 (Case 2)
    count_non_admissible(2, 2)

    # Example 4: a=3, b=4 (Case 2)
    count_non_admissible(3, 4)

if __name__ == "__main__":
    main()