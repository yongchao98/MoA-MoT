import math

def print_loss_formula():
  """
  This function prints the normalized AC loss formula for a superconductor
  in the critical state with an elliptical cross-section.
  """
  
  # The normalized variable i is the ratio of current amplitude to critical current
  i = "i"
  
  # The formula is derived from Norris's work for transport current losses.
  # Q = (mu_0 * Ic^2 / pi) * [ (1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2 ]
  # We need to calculate 2*pi*Q / (mu_0 * Ic^2)
  #
  # Normalized Loss = (2 * math.pi / (mu_0 * Ic**2)) * Q
  # Normalized Loss = (2 * math.pi) * (1 / math.pi) * [ ... ]
  # Normalized Loss = 2 * [ ... ]
  
  # Building the final string to print
  equation = f"2 * ((1 - {i}) * ln(1 - {i}) + (1 + {i}) * ln(1 + {i}) - {i}**2)"
  
  print("The normalized AC loss per cycle per unit length, 2*pi*Q / (mu_0 * Ic^2),")
  print("for an elliptic superconductor with i = Im/Ic < 1 is:")
  print("\nNormalized Loss = " + equation + "\n")
  print("Where:")
  print("  - 'i' is the normalized current amplitude (Im/Ic).")
  print("  - 'ln' is the natural logarithm.")

# Execute the function to print the formula
print_loss_formula()
