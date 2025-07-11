def get_stabilizing_controller_form():
  """
  This function prints the general form of all proper stabilizing controllers H_2(s)
  for the plant H_1(s) = s / (s^2 - 1).

  The controller is parametrized by a function K(s), which represents any stable,
  proper rational function satisfying the condition K(inf) != 1.
  """

  # Define the numerator and denominator as strings
  # H_2(s) = (X + D*K) / (Y - N*K)
  # Numerator = X + D*K = 4 + ((s-1)/(s+1)) * K
  # Denominator = Y - N*K = (s-1)/(s+1) - (s/(s+1)^2) * K
  # After simplification (see derivation in thought process):
  # Numerator(s) = (s+1) * [4*(s+1) + (s-1)*K(s)]
  #              = (s+1) * [4*s + 4 + (s-1)*K(s)]
  # Denominator(s) = (s-1)(s+1) - s*K(s)
  #                = s^2 - 1 - s*K(s)

  num_expression = "(s + 1) * ( (1*s - 1)*K(s) + 4*s + 4 )"
  den_expression = "1*s**2 - 1 - 1*s*K(s)"

  print("The set of all proper stabilizing controllers H_2(s) is given by:")
  print("H_2(s) = Num(s) / Den(s)\n")
  print("where K(s) is any stable and proper rational function such that its value at infinity, K(inf), is not equal to 1.\n")
  print(f"Num(s) = {num_expression}")
  print(f"Den(s) = {den_expression}")

get_stabilizing_controller_form()