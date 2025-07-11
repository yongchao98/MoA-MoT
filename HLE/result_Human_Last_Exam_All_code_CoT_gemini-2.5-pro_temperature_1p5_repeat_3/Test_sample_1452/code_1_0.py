def find_born_equation():
  """
  This function provides the famous equation associated with Max Born.

  Max Born's most significant contribution to quantum mechanics is the 'Born rule',
  which provides a physical interpretation for the wave function (Ψ). The rule states
  that the probability density of finding a particle at a specific location is
  proportional to the square of the magnitude of the wave function at that point.
  """
  name = "Max Born"
  # The equation is |Ψ|², where the '2' represents the square.
  equation_latex = r"\(|\Psi|^2\)"

  print(f'Name: "{name}"')
  print(f"Equation: {equation_latex}")

find_born_equation()