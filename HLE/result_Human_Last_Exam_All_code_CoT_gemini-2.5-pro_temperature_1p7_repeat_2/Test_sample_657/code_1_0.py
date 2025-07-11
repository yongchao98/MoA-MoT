def solve_controller():
  """
  This function formats and prints the general expression for the set of 
  all stabilizing controllers H_2(s).
  The expression is derived using the Youla-Kucera parametrization.
  """

  # Coefficients from the particular solution Yp(s) = 4s + 4
  y1 = 4
  y0 = 4

  # Coefficients from the particular solution Xp(s) = s - 1
  x1 = 1
  x0 = -1

  # Denominator of H1(s) is s^2 - 1
  d2 = 1
  d0 = -1
  
  # Print the numerator of H_2(s)
  print("The set of all stabilizing controllers H_2(s) is given by:")
  print("H_2(s) = N(s)/D(s), where:")
  
  # Format and print the numerator and denominator polynomials
  # N(s) = Y(s) = Y_p(s) - K(s)*D1(s)
  # D(s) = X(s) = X_p(s) + K(s)*N1(s)
  # Remember H1(s) = N1(s)/D1(s) = s / (s^2-1)
  
  print(f"N(s) = ({y1})s + ({y0}) - K(s)(s^2 + ({d0}))")
  print(f"D(s) = ({x1})s + ({x0}) + K(s)s")
  print("\nHere, K(s) is a parameter that represents any stable, proper rational function ensuring the controller H_2(s) itself is proper.")
  print("If K(s) is treated as an arbitrary polynomial as per the prompt, properness of H_2(s) is generally not guaranteed.")


solve_controller()
