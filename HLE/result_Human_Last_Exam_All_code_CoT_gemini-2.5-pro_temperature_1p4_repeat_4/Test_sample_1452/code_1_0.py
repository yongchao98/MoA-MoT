import re

def print_scientist_equation():
  """
  This function prints the name and famous equation for Max Born.
  The Born rule is a fundamental postulate of quantum mechanics that gives the
  probability that a measurement of a quantum system will yield a given result.
  """
  name = "Max Born"
  
  # The Born rule equation in LaTeX format.
  # P(x) is the probability density of finding the particle at position x.
  # |ψ(x)|² is the square of the absolute value of the wave function ψ at position x.
  equation_latex = r"P(x) = |\psi(x)|^2"

  # Extract numbers from the equation string to fulfill the requirement.
  # "Remember in the final code you still need to output each number in the final equation!"
  # The number in this equation is '2'. We will print it as part of the equation itself.
  
  print(f'Name: "{name}"')
  print(f"Equation: \\( {equation_latex} \\)")

print_scientist_equation()