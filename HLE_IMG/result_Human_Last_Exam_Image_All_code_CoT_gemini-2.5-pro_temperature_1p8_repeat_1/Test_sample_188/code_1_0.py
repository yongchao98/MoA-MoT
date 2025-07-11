def solve_chromatic_polynomial():
  """
  This function prints the calculated chromatic polynomial for the given graph.
  """
  # The chromatic polynomial P(k) was calculated using the deletion-contraction method.
  # The factored form is: P(k) = k * (k - 1) * (k - 2) * (k^2 - 4k + 5)
  # The expanded form is:
  
  # Coefficients for the polynomial P(k) = c5*k^5 + c4*k^4 + c3*k^3 + c2*k^2 + c1*k
  c5 = 1
  c4 = -7
  c3 = 19
  c2 = -23
  c1 = 10
  
  # Print the equation term by term as requested.
  # The 'k^5' term has a coefficient of 1, which is usually omitted in standard notation.
  print("The chromatic polynomial for the given graph is P(k), where:")
  print(f"P(k) = k^5 - {abs(c4)}*k^4 + {c3}*k^3 - {abs(c2)}*k^2 + {c1}*k")

solve_chromatic_polynomial()