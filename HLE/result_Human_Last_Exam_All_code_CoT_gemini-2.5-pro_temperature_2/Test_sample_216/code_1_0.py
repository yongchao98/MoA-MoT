def print_performance_bound_formula():
  """
  This function prints the symbolic formula for the tightest upper bound
  on the performance difference J(pi^*) - J(pi_hat).
  """
  # Symbolic variables in the formula
  H = "H" # Horizon
  R_max = "R_max" # Maximum reward
  abs_A = "|A|" # Size of the action space
  lambda_param = "lambda" # Hyperparameter

  # The derived bound formula includes the horizon squared, reflecting the compounding
  # error in imitation learning. The bound on the total variation (TV) risk is
  # directly substituted into the standard performance difference bound.
  # J(pi^*) - J(pi_hat) <= H^2 * R_max * TV_risk
  # J(pi^*) - J(pi_hat) <= H**2 * R_max * |A| * (1 - exp(-lambda))
  
  # Printing each term clearly, as requested.
  term1 = f"{H}**2"
  term2 = R_max
  term3 = abs_A
  term4 = f"(1 - exp(-{lambda_param}))"

  print(f"The tightest upper bound for J(pi^*) - J(pi_hat) is: {term1} * {term2} * {term3} * {term4}")

print_performance_bound_formula()