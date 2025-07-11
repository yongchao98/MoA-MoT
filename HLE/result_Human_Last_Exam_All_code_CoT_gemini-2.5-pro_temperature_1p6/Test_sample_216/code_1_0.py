import math

# The parameters of the problem are symbolic:
# H: horizon
# |A|: size of the action space
# lambda: a hyperparameter

# The tightest upper bound on the performance difference J(pi^*) - J(pi_hat)
# is derived from imitation learning theory. The final formula is a product of
# the horizon H and the given upper bound on the population total variation risk.

# We will print the formula symbolically.
H = "H"
A_size = "|A|"
hyperparam_lambda = "lambda"

# The numbers in the equation are 1 and -1.
constant_one = 1
constant_minus_one = -1

# Print the final expression for the tightest upper bound.
print(f"The tightest upper bound for J(pi^*) - J(pi_hat) is:")
print(f"{H} * {A_size} * ({constant_one} - exp({constant_minus_one} * {hyperparam_lambda}))")