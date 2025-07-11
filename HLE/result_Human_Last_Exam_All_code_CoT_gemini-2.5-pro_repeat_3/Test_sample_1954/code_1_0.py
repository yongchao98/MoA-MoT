import math

# This script provides the formula for the minimax risk based on the derivation.

# Problem interpretation:
# We have n i.i.d. observations X_i ~ Bin(n, theta).
# The sufficient statistic is S = sum(X_i), which follows Bin(N, theta) where N = n*n.
# The loss function is squared error: L(d, theta) = (d - theta)^2.

# The minimax risk for estimating theta from a Bin(N, theta) observation
# is found using a Bayes estimator with constant risk. This is achieved
# by using a Beta(alpha, beta) prior with alpha = beta = sqrt(N)/2.

# The resulting minimax risk is given by the formula:
# R_minimax = 1 / (4 * (sqrt(N) + 1)^2)

# Substituting N = n^2, the formula becomes:
# R_minimax = 1 / (4 * (sqrt(n^2) + 1)^2) = 1 / (4 * (n + 1)^2)

# The following code prints the final formula, explicitly showing each number.

numerator = 1
denominator_coeff = 4
inner_term_constant = 1

print("The minimax risk for estimating theta is given by the formula:")
print(f"{numerator} / ({denominator_coeff} * (n + {inner_term_constant})^2)")
