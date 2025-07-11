import math
from scipy.stats import norm

# Given parameters
gamma = 2.1  # Power-law exponent
alpha = 2.5  # Pareto distribution shape
epsilon = 0.05  # Marginal completeness tolerance
confidence = 0.99  # Confidence level

# Calculate delta (significance level)
delta = 1 - confidence

# Calculate the Z-score for the two-sided confidence interval
# We use 1 - delta/2 because the confidence interval is two-sided
Z = norm.ppf(1 - delta / 2)

# Calculate the heterogeneity factor based on gamma and alpha
heterogeneity_factor = (gamma - 2) * (alpha - 2)

# Calculate the numerator of the formula, based on Chernoff-Hoeffding bounds
numerator = math.log(2 / delta)

# Calculate the denominator of the formula
# This combines the statistical requirement (Z/epsilon)^2 and the network complexity
denominator = ((Z / epsilon)**2) * heterogeneity_factor

# Calculate the minimum sampling ratio r
r = numerator / denominator

# Output the equation and the final result
print("The formula for the sampling ratio 'r' is derived by combining statistical requirements with network complexity factors:")
print("r = log(2 / delta) / ( (Z / epsilon)^2 * (gamma - 2) * (alpha - 2) )")
print("\nWhere:")
print(f"gamma = {gamma}")
print(f"alpha = {alpha}")
print(f"epsilon = {epsilon}")
print(f"confidence = {confidence} -> delta = {delta:.2f}")
print(f"Z-score for {confidence*100}% confidence = {Z:.4f}")
print("\nCalculation:")
print(f"r = log(2 / {delta:.2f}) / ( ({Z:.4f} / {epsilon})**2 * ({gamma} - 2) * ({alpha} - 2) )")
print(f"r = {numerator:.4f} / ( ({Z/epsilon:.4f})**2 * {gamma-2:.1f} * {alpha-2:.1f} )")
print(f"r = {numerator:.4f} / ( { (Z/epsilon)**2 :.4f} * {heterogeneity_factor:.2f} )")
print(f"r = {numerator:.4f} / {denominator:.4f}")
print(f"r = {r:.4f}")
print("\nTherefore, the minimum required ratio r is:")
print(f"{r:.4f}")
<<<0.0399>>>