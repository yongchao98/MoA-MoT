import math

def max_complex_zeros(N):
  """
  Calculates the maximal possible number of complex zeros for the given matrix B(k).
  The formula is derived by finding the degree of the final polynomial equation for k_1.

  The derivation steps:
  1. The determinant equation det(B(k))=0 is polynomial in k_1 and radicals sqrt(k_1^2 + Delta_j).
  2. A systematic elimination of the radicals leads to a single polynomial in k_1.
  3. Let y = k_1. The degree of the polynomial in y is tracked through the elimination process.
  4. Let d_k be the maximal degree in y of the coefficients of the equation written as a polynomial in s_k = sqrt(y^2 + Delta_k).
  5. A recurrence relation for the degree is found: d_{k-1} = 2*d_k + 4.
  6. The base case for the recurrence is d_{N-1} = 4.
  7. The degree of the final polynomial in y is found to be 2^{N+1} - 4.
  8. This degree is the maximal number of roots.
  """
  if not isinstance(N, int) or N < 2:
    return "N must be an integer greater than or equal to 2."
  
  # The final equation for the number of zeros is 2^(N+1) - 4
  result = 2**(N + 1) - 4
  return result

# We want to print the general formula as the answer.
# The question also mentions outputting each number in the final equation.
# The formula is 2^(N+1) - 4
# Let's show this for N.
N_str = "N"
final_equation = f"2^({N_str}+1) - 4"
print(f"The general formula for the maximal number of complex zeros is: {final_equation}")

# Example calculations for specific N values
print("\nFor specific values of N:")
for n_val in range(2, 6):
    num_zeros = max_complex_zeros(n_val)
    # The final equation for a given N value.
    equation = f"2^({n_val}+1) - 4 = {2**(n_val+1)} - 4 = {num_zeros}"
    print(f"For N = {n_val}, the maximal number of zeros is given by the equation: {equation}")