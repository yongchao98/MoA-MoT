import math

def calculate_max_complex_zeros(N):
  """
  Calculates the maximal possible number of complex zeros for the matrix B(k).

  The number of zeros is given by the formula N * 2^(N-1), derived from
  analyzing the system of polynomial equations that define the problem.
  By choosing the matrix A appropriately, all these zeros can be made "complex"
  in the sense required by the problem (Re(k_j) != 0 and Im(k_j) != 0).

  This function demonstrates the calculation for a given N.
  """
  if not isinstance(N, int) or N <= 0:
    print("Error: N must be a positive integer.")
    return

  # The formula for the maximal number of zeros is N * 2^(N-1)
  max_zeros = N * (2**(N-1))

  # Print the result showing the final equation with numbers.
  print(f"For a matrix of size N = {N}:")
  print(f"The maximal number of complex zeros is calculated as {N} * 2^({N}-1) = {max_zeros}")

# Example calculations for N = 2, 3, and 4.
calculate_max_complex_zeros(2)
print("-" * 20)
calculate_max_complex_zeros(3)
print("-" * 20)
calculate_max_complex_zeros(4)
