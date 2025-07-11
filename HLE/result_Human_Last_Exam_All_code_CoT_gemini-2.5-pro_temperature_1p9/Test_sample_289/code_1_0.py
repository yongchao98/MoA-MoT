import cmath
import numpy as np

def find_largest_set_size():
    """
    Finds the largest size |S| of a set S of non-real eigenvalues
    for a matrix A satisfying A^3 = A^*.

    The core of the problem is to solve the equation lambda^3 = conjugate(lambda),
    which is the condition that any eigenvalue lambda must satisfy.
    We solve this by expressing lambda in polar coordinates: lambda = r * exp(i*theta).
    The equation becomes r^3 * exp(i*3*theta) = r * exp(-i*theta).
    This leads to two conditions:
    1. For the magnitude r: r^3 = r  => r=0 or r=1.
    2. For the angle theta (when r=1): 3*theta = -theta + 2*k*pi => theta = k*pi/2 for integer k.
    """

    possible_eigenvalues = set()

    # Case 1: r=0
    # This gives the trivial solution lambda = 0, which is real.
    lambda_zero = 0.0 + 0.0j
    possible_eigenvalues.add(lambda_zero)

    # Case 2: r=1
    # We find all unique solutions for theta = k*pi/2 by checking k=0, 1, 2, 3.
    # Higher integer values of k will only produce repeated solutions.
    for k in range(4):
        theta = k * cmath.pi / 2
        lambda_k = cmath.exp(1j * theta)
        # Round the result to handle potential floating point inaccuracies.
        # This helps ensure that complex numbers like 1.0 + 0.0j are correctly identified.
        rounded_lambda = round(lambda_k.real, 10) + round(lambda_k.imag, 10) * 1j
        possible_eigenvalues.add(rounded_lambda)

    # The set S consists of non-real eigenvalues. A number is non-real if its imaginary part is non-zero.
    S = {val for val in possible_eigenvalues if val.imag != 0}

    # The largest possible size of a set of non-real eigenvalues is the total count
    # of all possible non-real eigenvalues. We can construct a diagonal matrix with these
    # eigenvalues which satisfies the condition.
    max_size = len(S)

    print("Let the size of the set be |S|.")
    print("The final equation for the size of the set is: |S| =", max_size)
    print("The numbers that are elements of this largest possible set S are:")
    for val in sorted(list(S), key=lambda x: x.imag):
        print(val)
    
    print("\nThe largest possible size |S| is:")
    print(max_size)

if __name__ == "__main__":
    find_largest_set_size()