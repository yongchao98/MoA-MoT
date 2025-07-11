# Initial number of good (white) and defective (black) products
W0 = 2
B0 = 1

# The problem is to find the probability p that the number of good and defective products
# eventually become equal.
# This probability can be determined using the properties of the Beta distribution,
# which describes the limit of the proportion of good products.

# The limiting proportion Y_inf follows a Beta(W0, B0) distribution.
# For our case, Y_inf ~ Beta(2, 1).
# The PDF of Beta(2, 1) is f(y) = 2y for y in [0, 1].

# Let B be the event that the limiting proportion Y_inf is >= 1/2.
# P(B) = integral from 1/2 to 1 of 2y dy = [y^2] from 1/2 to 1.
P_B = 1**2 - (1/2)**2

# The probability P(B) can also be expressed in terms of p:
# P(B) = P(Y_inf >= 1/2 | T < inf)*p + P(Y_inf >= 1/2 | T = inf)*(1-p)
# P(Y_inf >= 1/2 | T < inf) = 1/2 (due to symmetry of the subsequent process)
# P(Y_inf >= 1/2 | T = inf) = 1 (since T=inf implies Y_t > 1/2 for all t)
# So, P(B) = (1/2)*p + 1*(1-p) = 1 - p/2.

# Now we can solve for p:
# P_B = 1 - p/2
# p/2 = 1 - P_B
# p = 2 * (1 - P_B)
prob_p = 2 * (1 - P_B)

# The result is the exact probability, which is also the least upper bound.
print(f"The initial number of good products is W_0 = {W0}")
print(f"The initial number of defective products is B_0 = {B0}")
print("The limiting proportion of good products follows a Beta(2, 1) distribution.")
print("Let B be the event that this limiting proportion is at least 1/2.")
print(f"The probability of event B is P(B) = 1^2 - (1/2)^2 = {P_B}")
print("The probability p is related to P(B) by the equation: P(B) = 1 - p/2")
print(f"Solving for p: p = 2 * (1 - P(B))")
print(f"p = 2 * (1 - {P_B}) = {prob_p}")
print(f"The upper bound for the probability is {prob_p}")
