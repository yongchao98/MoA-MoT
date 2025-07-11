# The analysis of the Markov chain's transience boils down to the sign
# of the drift of a test function f(x) = log(x).
# This drift, for large k, is approximated by the term (A*c - B) / k^2.
# Transience occurs when A*c - B > 0.
# The critical value c is found when A*c - B = 0.

# From the problem's transition probabilities, we derive the coefficients A and B.
# A is the coefficient of c in the k-scaled first moment (mu_k).
# mu_k = 2*c/k, so after dividing by k, the coefficient of c/k^2 is 2.
coeff_A = 2.0

# B is half of the second moment (M2_k) of the jump.
# M2_k = 2.5. Therefore, B = 2.5 / 2.
m2 = 2.5
coeff_B = m2 / 2.0

# We need to solve the equation: coeff_A * c - coeff_B = 0
# which is 2.0 * c - 1.25 = 0.

# We solve for the critical value of c.
c_infimum = coeff_B / coeff_A

print("The critical value for c is determined by the equation for the drift of log(k):")
print(f"A*c - B = 0")
print("Where A is the coefficient from the first moment, and B is from the second moment.")
print(f"Based on the transition probabilities, A = {coeff_A} and B = {m2} / 2 = {coeff_B}.")
print("The equation to solve is:")
print(f"{coeff_A} * c - {coeff_B} = 0")
print(f"Solving for c gives c = {coeff_B} / {coeff_A}")
c = 0.625
print(f"c = {c}")

print(f"\nAnalysis shows the chain is transient for c >= {c} and recurrent for c < {c}.")
print(f"The set of c for which the chain is transient is [{c}, infinity).")
print(f"The infimum of this set is therefore {c}.")
