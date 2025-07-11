import math

def final_equation():
  """
  This function prints the derived formula for the instantaneous force f_x(t) on the coil.
  The final expression corresponds to choice D.
  """
  # The formula is constructed as a string. The numbers 2 and 1 are explicitly included
  # as requested.
  equation_string = "f_x(t) = 2 * pi * R * N * (mu_0 * (1 - alpha_T * (T - T_0)) * N_0 * I_0 * i_0 * sin(omega * t)) / (g**2 * (1 + (mu_0 * N_0 * I_0) / (g * B_s)))"
  
  print("The final equation for the instantaneous force is:")
  print(equation_string)

final_equation()