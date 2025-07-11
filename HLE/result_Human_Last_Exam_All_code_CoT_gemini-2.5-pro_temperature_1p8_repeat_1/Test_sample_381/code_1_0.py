import math

def get_upper_bound_expression():
  """
  Calculates and prints the expression for the upper bound.
  """
  # The upper bound for ||B * Q_{0,M}||_infinity is derived as follows:
  # 1. Using the submultiplicative property of matrix norms:
  #    ||B * Q||_inf <= ||B||_inf * ||Q||_inf

  # 2. Bounding ||Q||_inf:
  #    Q is a product of matrices of the form DP.
  #    The infinity norm of DP is max_i(sum_j |(DP)_ij|) = max_i(D_ii * sum_j(P_ij)) = max_i(D_ii * 1) <= 1.
  #    Since ||DP||_inf <= 1, the norm of the product Q is also ||Q||_inf <= 1.

  # 3. Bounding ||B||_inf:
  #    B is an (N-1)xN matrix with orthonormal rows b_i, where ||b_i||_2 = 1.
  #    By the Cauchy-Schwarz inequality, a row's l1-norm ||b_i||_1 <= sqrt(N) * ||b_i||_2 = sqrt(N).
  #    Thus, ||B||_inf = max_i(||b_i||_1) <= sqrt(N).

  # 4. Combining the bounds:
  #    ||B * Q||_inf <= ||B||_inf * ||Q||_inf <= sqrt(N) * 1 = sqrt(N).

  # The equation for the upper bound is 1 * sqrt(N).
  # The number in the final equation is 1.
  number_in_equation = 1
  
  # The variable N is symbolic, so we print the expression as a string.
  print(f"The upper-bound is: {number_in_equation} * sqrt(N)")

get_upper_bound_expression()