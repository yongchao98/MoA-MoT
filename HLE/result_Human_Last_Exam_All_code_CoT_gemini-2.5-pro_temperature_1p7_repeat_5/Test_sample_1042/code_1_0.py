import sys

# Suppose for a given value L=l, we have the following quantities identified from data.
# This corresponds to the regression of Y on A at a specific value of L.
E_Y_given_A1_L = 10.0  # E(Y | A=1, L=l)
E_Y_given_A0_L = 5.0   # E(Y | A=0, L=l)

# This corresponds to the propensity score at L=l.
p_A1_given_L = 0.6     # P(A=1 | L=l)
p_A0_given_L = 1.0 - p_A1_given_L

# These are the average counterfactual outcomes, which the problem states we can identify.
# We choose values that are consistent with the other quantities.
E_Y1_given_L = 8.5 # E(Y^1 | L=l)
E_Y0_given_L = 6.5 # E(Y^0 | L=l)

# --- Identification Step-by-Step ---

print("Goal: Identify all conditional expectations E(Y^a | A=a', L=l)\n")

# Case 1: a = a' (Identified by consistency)
# Calculate E(Y^1 | A=1, L=l)
E_Y1_given_A1_L = E_Y_given_A1_L
print("1. Identifying E(Y^a=1 | A=1, L=l):")
print("   By consistency, E(Y^1 | A=1, L=l) = E(Y | A=1, L=l)")
print(f"   E(Y^1 | A=1, L=l) = {E_Y1_given_A1_L}\n")

# Calculate E(Y^0 | A=0, L=l)
E_Y0_given_A0_L = E_Y_given_A0_L
print("2. Identifying E(Y^a=0 | A=0, L=l):")
print("   By consistency, E(Y^0 | A=0, L=l) = E(Y | A=0, L=l)")
print(f"   E(Y^0 | A=0, L=l) = {E_Y0_given_A0_L}\n")

# Case 2: a != a' (Identified by Law of Total Expectation)
# Check for positivity violation before proceeding
if p_A0_given_L == 0 or p_A1_given_L == 0:
    print("Positivity violation: Cannot identify due to zero probability.", file=sys.stderr)
    sys.exit(1)

# Calculate E(Y^1 | A=0, L=l)
print("3. Identifying E(Y^a=1 | A=0, L=l):")
print("   From the law of total expectation: E(Y^1|L) = E(Y^1|A=1,L)P(A=1|L) + E(Y^1|A=0,L)P(A=0|L)")
print("   Solving for E(Y^1|A=0,L) gives: [E(Y^1|L) - E(Y^1|A=1,L)P(A=1|L)] / P(A=0|L)")
print("   Substituting identifiable quantities:")
# Using a temp variable for clarity
numerator = E_Y1_given_L - E_Y1_given_A1_L * p_A1_given_L
E_Y1_given_A0_L = numerator / p_A0_given_L
print(f"   E(Y^1|A=0,L) = [{E_Y1_given_L} - {E_Y1_given_A1_L} * {p_A1_given_L}] / {p_A0_given_L}")
print(f"   E(Y^1 | A=0, L=l) = {E_Y1_given_A0_L}\n")


# Calculate E(Y^0 | A=1, L=l)
print("4. Identifying E(Y^a=0 | A=1, L=l):")
print("   From the law of total expectation: E(Y^0|L) = E(Y^0|A=1,L)P(A=1|L) + E(Y^0|A=0,L)P(A=0|L)")
print("   Solving for E(Y^0|A=1,L) gives: [E(Y^0|L) - E(Y^0|A=0,L)P(A=0|L)] / P(A=1|L)")
print("   Substituting identifiable quantities:")
numerator = E_Y0_given_L - E_Y0_given_A0_L * p_A0_given_L
E_Y0_given_A1_L = numerator / p_A1_given_L
print(f"   E(Y^0|A=1,L) = [{E_Y0_given_L} - {E_Y0_given_A0_L} * {p_A0_given_L}] / {p_A1_given_L}")
print(f"   E(Y^0 | A=1, L=l) = {E_Y0_given_A1_L}\n")