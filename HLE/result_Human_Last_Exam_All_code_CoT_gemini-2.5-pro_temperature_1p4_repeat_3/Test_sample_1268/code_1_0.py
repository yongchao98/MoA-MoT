import math

def is_squarefree(n):
    """
    Checks if a positive integer is squarefree.
    A number is squarefree if it is not divisible by any perfect square other than 1.
    """
    if n <= 0:
        return False
    i = 2
    # Check for divisibility by squares of primes up to sqrt(n)
    while i * i <= n:
        if n % (i * i) == 0:
            return False
        i += 1
    return True

def calculate_norm_bound():
    """
    Calculates and prints the upper bound for the maximum norm for a given
    squarefree natural number N.
    """
    try:
        n_input = input("Please enter a squarefree natural number N: ")
        N = int(n_input)
        if N <= 0:
            raise ValueError("N must be a positive integer.")
    except (ValueError, TypeError):
        print(f"Invalid input. Please enter a positive integer.")
        return

    if not is_squarefree(N):
        print(f"Error: The number {N} is not squarefree. Please provide a number that is not divisible by any perfect square other than 1.")
        return

    # For a real quadratic field Q(sqrt(N)), the dimension of the lattice is n=2.
    n = 2

    # The discriminant Delta_K of the field Q(sqrt(N)) depends on N mod 4.
    if N % 4 == 1:
        delta_K = N
    else:  # This covers N % 4 == 2 and N % 4 == 3
        delta_K = 4 * N
    
    # The covolume V of the lattice is the square root of the absolute value of the discriminant.
    V = math.sqrt(delta_K)
    
    # According to Minkowski's theorem, an upper bound for the infinity norm
    # of the shortest non-zero vector is V^(1/n).
    upper_bound = V**(1/n)
    
    # Print the explanation and the results
    print(f"\nFor the number field associated with the squarefree number N = {N}:")
    print(f"1. The dimension of the lattice is n = {n}.")
    print(f"2. The field discriminant is Delta_K = {delta_K}.")
    print(f"3. The lattice covolume is V = sqrt({delta_K}) = {V:.5f}.")
    
    print("\nAn upper bound for the maximum norm (k_k,inf) is given by V^(1/n).")

    # Final equation with numbers, as requested.
    print("\nThe final equation for the upper bound is:")
    print(f"k_{{k,inf}} <= {V:.5f} ^ (1 / {n})")
    print(f"k_{{k,inf}} <= {upper_bound:.5f}")


if __name__ == '__main__':
    calculate_norm_bound()