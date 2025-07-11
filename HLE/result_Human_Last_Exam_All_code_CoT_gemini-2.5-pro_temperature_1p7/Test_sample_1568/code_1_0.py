import math

def display_product_formula():
  """
  Prints the closed-form expression for the given infinite product.
  Note: This code provides the symbolic answer, not a numerical computation.
  The symbols 'Gamma', 'exp', 'pi', and 'i' represent their usual mathematical meanings.
  """
  
  # The product starts from n=3, so we divide the general formula
  # by the terms for n=1 and n=2.
  n1_term = "1 - z^3/1^3"
  n2_term = "1 - z^3/2^3"
  
  # The product of Gamma functions in the denominator of the general formula.
  # Using w = exp(2*pi*i/3) notation is also common.
  gamma_terms = "Gamma(1 - z) * Gamma(1 - exp(2*pi*i/3)*z) * Gamma(1 - exp(4*pi*i/3)*z)"

  # Constructing the final expression string.
  # The numbers from the original equation (3) and the derivation (1, 2) are present.
  final_expression = f"1 / ( ({n1_term}) * ({n2_term}) * ({gamma_terms}) )"
  
  print(final_expression)

display_product_formula()