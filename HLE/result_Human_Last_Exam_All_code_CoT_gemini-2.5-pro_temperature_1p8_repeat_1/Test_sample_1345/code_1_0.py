import sys

def solve():
    """
    This function calculates the maximal possible number of complex zeros for the
    given matrix equation.

    The problem is to find the maximum number of complex solutions k_j to the
    equation det(B(k)) = 0, where k_j are related by k_j^2 = k_1^2 + Delta_j.
    A solution is "complex" if k_j is neither real nor purely imaginary.

    This problem can be transformed into finding the non-real roots of a
    polynomial P(z) in the variable z = k_1^2. The degree of this polynomial
    is D = N * 2^(N-1).

    A solution k_j is complex if and only if z = k_1^2 is non-real.
    Since P(z) has real coefficients, its non-real roots come in conjugate pairs.

    - If N=1, the degree D = 1*2^0 = 1. A real polynomial of odd degree must
      have at least one real root. The maximum number of non-real roots is D-1=0.

    - If N>=2, the degree D = N * 2^(N-1) is always even. It is possible
      to choose the matrix A such that all D roots are non-real.

    Therefore, the maximal number of complex zeros is 0 for N=1 and N * 2^(N-1)
    for N >= 2.
    """
    
    # Prompting the user to enter the value for N
    try:
        n_input = input("Enter the integer value for N (the size of the matrix): ")
        N = int(n_input)
        if N <= 0:
            print("N must be a positive integer.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer for N.")
        return
    except Exception as e:
        print(f"An error occurred: {e}")
        return

    if N == 1:
        max_zeros = 0
    else:
        # Using integer arithmetic for precision with large numbers
        power_of_2 = 1 << (N - 1) # Equivalent to 2**(N-1)
        max_zeros = N * power_of_2

    print(f"For N = {N}:")
    print(f"The degree of the polynomial in z=k_1^2 is D = N * 2^(N-1).")
    
    if N == 1:
        print(f"D = 1 * 2^(1-1) = 1. The degree is odd, so there must be at least one real root.")
        print(f"The maximal number of complex zeros is 0.")
    else:
        print(f"D = {N} * 2^({N-1}) = {N} * {1 << (N-1)} = {max_zeros}.")
        print(f"Since D is even, it's possible for all roots to be non-real.")
        print(f"The maximal number of complex zeros is {max_zeros}.")
        
    print("\n<<<{}>>>".format(max_zeros))


solve()