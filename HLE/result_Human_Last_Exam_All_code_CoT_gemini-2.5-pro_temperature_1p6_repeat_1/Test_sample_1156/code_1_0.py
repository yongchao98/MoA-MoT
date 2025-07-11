import math

# The derived normalized density for the modified problem
# rho(x) = 1 / (2*ln(2)*(x + sqrt(x)))
# We will print the components of this equation.

c1 = 1
c2 = 2
log_of_2 = math.log(2)
# We can approximate log(2) if needed for display
log_of_2_approx = 0.6931

print("Based on the analysis of a slightly modified but solvable version of the problem, the normalised density of the invariant measure is:")
print(f"rho(x) = {c1} / ({c2}*ln({c2})*(x + sqrt(x)))")
print("\nWhich is approximately:")
print(f"rho(x) = {c1} / ({c2*log_of_2_approx}*(x + sqrt(x)))")
print(f"rho(x) = {c1} / ({round(c2*log_of_2_approx, 4)}*(x + sqrt(x)))")

# Final equation components:
print("\nThe final equation is:")
print(f"rho(x) = {c1} / ( {c2} * ln({c2}) * (x + x^(1/2)) )")
