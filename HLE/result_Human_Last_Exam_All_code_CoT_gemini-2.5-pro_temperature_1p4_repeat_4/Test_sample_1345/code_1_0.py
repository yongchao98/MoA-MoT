import sys

def max_complex_zeros(N):
  """
  Calculates the maximal possible number of complex zeros based on the derived formula.
  A zero k_j is complex if Re(k_j) != 0 and Im(k_j) != 0.
  
  Args:
    N: The dimension of the matrix.
    
  Returns:
    The maximal number of complex zeros.
  """
  if not isinstance(N, int) or N <= 0:
    print("Error: N must be a positive integer.", file=sys.stderr)
    return None
  
  if N == 1:
    # For N=1, the degree of the characteristic polynomial for k_1 is 1.
    # A polynomial of odd degree with real coefficients must have at least one real root.
    # This real root for k_1 does not satisfy the condition for a complex zero.
    # Therefore, the maximum number of complex zeros is 0.
    return 0
  else:
    # For N > 1, the degree of the polynomial for k_1 is d = N * 2^(N-1).
    # Since N-1 >= 1, 2^(N-1) is even, making d even.
    # For a polynomial of even degree, it is possible to choose the matrix A
    # such that all roots are complex and not purely imaginary.
    # The maximal number of such zeros is equal to the degree of the polynomial.
    return N * (2**(N-1))

if __name__ == "__main__":
    print("The maximal number of complex zeros is given by the formula:")
    print("f(N) = 0, if N = 1")
    print("f(N) = N * 2^(N-1), if N > 1")
    print("\nCalculating for specific values of N:")
    
    for n_val in range(1, 6):
        power = n_val - 1
        base = 2
        result = max_complex_zeros(n_val)
        
        if n_val == 1:
            print(f"For N = {n_val}: The result is 0.")
        else:
            # Output each number in the final equation
            print(f"For N = {n_val}: The result is {result}. Calculated as: {n_val} * {base}^({power}) = {result}")
