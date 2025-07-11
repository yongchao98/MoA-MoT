def print_loss_formula():
  """
  This function prints the normalized AC loss formula for an elliptical superconductor.
  """
  # The formula for the normalized AC loss 2*pi*Q/(mu_0 * Ic^2)
  # as a function of i = Im/Ic, for i < 1.
  # 'ln' denotes the natural logarithm.
  # '**2' denotes squaring.
  formula = "2 * [ (1 - i)*ln(1 - i) + (1 + i)*ln(1 + i) - i**2 ]"
  
  print("The normalized AC loss per cycle per unit length, 2*pi*Q/(mu_0 * Ic^2), as a function of i = Im/Ic is given by the following expression:")
  print(formula)

# Execute the function to print the result.
print_loss_formula()