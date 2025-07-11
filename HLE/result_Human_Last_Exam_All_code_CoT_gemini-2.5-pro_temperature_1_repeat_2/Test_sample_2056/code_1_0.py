import sympy as sp

def get_final_expression():
  """
  This function provides the derived analytical expression for l_k(n).
  The derivation is explained in the text preceding this code block.
  """
  
  # Define the symbols n and k for the symbolic expression
  n, k = sp.symbols('n k')

  # Construct the expression based on the derivation
  # Term 1: from the determinant of Sigma
  term1 = sp.Rational(1, 2) * sp.log(n + 1)
  
  # Term 2: from the distance and Jacobian term
  term2 = (n - 1) * sp.log(k)
  
  # Term 3: from the quadratic form in the Gaussian exponent
  term3 = -k**2 * (2*n - 1) / n
  
  # The final expression for l_k(n)
  l_k_n = term1 + term2 + term3
  
  # The prompt asks to output each number in the final equation.
  # We print the complete formula as a string, which contains all constants and operations.
  final_equation_string = f"l_k(n) = {l_k_n}"

  print(final_equation_string)

if __name__ == '__main__':
  get_final_expression()