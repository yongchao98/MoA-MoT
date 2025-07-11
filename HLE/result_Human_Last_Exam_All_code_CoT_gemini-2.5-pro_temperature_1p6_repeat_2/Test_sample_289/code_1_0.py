import numpy as np

def solve_for_eigenvalues():
    """
    Solves for the largest size of a set of non-real eigenvalues for a matrix
    A where A^3 = A^*.
    """
    # The condition on the eigenvalues lambda is lambda^3 = conjugate(lambda).
    # In polar form (lambda = r * e^(i*theta)), this splits into two equations.

    # 1. Equation for the magnitude r: r^3 = r
    # This can be written as r^3 - r = 0.
    # The numbers in the polynomial r^3 - r are for the powers 3 and 1.
    print("Derived equation for magnitude r: r**3 = r")
    # The solutions for r >= 0 are r=0 and r=1.

    # 2. Equation for the argument theta: e^(i*3*theta) = e^(-i*theta)
    # This implies 4*theta = 2*k*pi for any integer k.
    # The numbers in this equation are 4 and 2.
    print("Derived equation for argument theta: 4 * theta = 2 * k * pi")
    print("-" * 30)

    # We now find all unique eigenvalues based on these equations.
    all_eigenvalues = set()

    # Case 1: r = 0
    # This gives the eigenvalue lambda = 0.
    all_eigenvalues.add(0.0)

    # Case 2: r = 1
    # We find the angles for k=0, 1, 2, 3 (other integers will repeat the solutions)
    # theta = k * pi / 2
    for k in range(4):
        theta = k * np.pi / 2
        # Calculate eigenvalue: lambda = 1 * e^(i*theta)
        eigenvalue = np.exp(1j * theta)
        # np.round is used to clean up floating point inaccuracies (e.g., e^(i*pi) = -1+1.22e-16j)
        all_eigenvalues.add(np.round(eigenvalue, 10))
    
    print(f"The set of all possible eigenvalues is: {all_eigenvalues}")

    # Now, we filter for the set S of non-real eigenvalues.
    # A complex number z is real if its imaginary part is zero.
    non_real_eigenvalues = {z for z in all_eigenvalues if z.imag != 0}

    print(f"The set S of non-real eigenvalues is: {non_real_eigenvalues}")

    # The size |S| is the number of elements in this set.
    size_S = len(non_real_eigenvalues)

    print(f"\nThe largest possible size |S| of the set of non-real eigenvalues is: {size_S}")

if __name__ == "__main__":
    solve_for_eigenvalues()
