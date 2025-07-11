# A simple demonstration of the rules of Berezin integration for a single Grassmann variable η.

# In this formalism, any function of η can be written as f(η) = A + Bη,
# because η² = 0 due to the anti-commuting nature of Grassmann variables,
# which encodes the Pauli exclusion principle.
# We can represent such a function by a tuple of its coefficients (A, B).

def berezin_integral(function_coeffs):
  """
  Calculates the Berezin integral for a function f(η) = A + Bη.
  The Berezin integral is defined by the rules ∫dη = 0 and ∫dη η = 1.
  Therefore, ∫dη (A + Bη) = A*∫dη + B*∫dη η = A*0 + B*1 = B.

  Args:
    function_coeffs: A tuple (A, B) representing the function coefficients.

  Returns:
    The result of the integral, which is the coefficient B.
  """
  A, B = function_coeffs
  return B

# --- Main Demonstration ---

print("The integration measure 'dη' for a Grassmann variable is defined by its action.")
print("Let's test the two fundamental rules of this integration.")
print("-" * 50)

# Rule 1: The integral of a constant is 0.
# The function is f(η) = 1, which means A=1, B=0.
constant_function = (1, 0)
integral_of_constant = berezin_integral(constant_function)
print("Demonstrating the first rule for f(η) = 1:")
print(f"The equation is: ∫ dη * 1 = {integral_of_constant}")
print("")


# Rule 2: The integral of the variable itself is 1.
# The function is f(η) = η, which means A=0, B=1.
linear_function = (0, 1)
integral_of_eta = berezin_integral(linear_function)
print("Demonstrating the second rule for f(η) = η:")
print("This rule provides the fundamental normalization value for the measure.")
print(f"The equation is: ∫ dη * η = {integral_of_eta}")
print("-" * 50)
print("This normalization value of 1 is the essential non-trivial value that defines the Grassmann integral measure.")
