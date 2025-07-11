import math

# Define the symbols used in the formula
# H: horizon of the episode
# abs_A: size of the discrete action space, |A|
# lambda_val: a hyperparameter of the algorithm

# We do not have numerical values for these, so we will represent them as strings.
H = "H"
abs_A = "|A|"
lambda_val = "lambda"

# The formula for the tightest upper bound of J(pi*) - J(pi_hat) is derived
# from standard imitation learning theory for behavioral cloning.
# The bound is H^2 * R_max * E[TV_distance], where R_max is assumed to be 1.
# The expectation of the TV distance is given by the TV risk, which is bounded by the expression provided.

# Construct the final formula string
# We represent H^2 as H**2 and e^-lambda as exp(-lambda) for clarity.
bound_expression = f"{H}**2 * {abs_A} * (1 - exp(-{lambda_val}))"

# Print the final result in the format of an equation
print(f"J(pi^*) - J(pi_hat) <= {bound_expression}")