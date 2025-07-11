# Step 1: Determine the cardinality of S/A.
# S/A is the quotient group *R/Z.
# The cardinality of the hyperreals *R is Beth_1.
# The cardinality of Z is Beth_0.
# The cardinality of the quotient |*R/Z| is Beth_1.
card_S_div_A = "Beth_1"

# Step 2: Determine the cardinality of B/S.
# B/S is the quotient group Z/Im(g), where g is the canonical map from *R to Z.
# The map g is the zero homomorphism. Its image is {0}.
# The cardinality of Z/{0} is |Z| = Beth_0.
card_B_div_S = "Beth_0"

# Step 3: Determine the cardinality of H_1(B/A, Q).
# The space B/A is Z/{0}, which is Z.
# With the discrete topology, H_1(Z, Q) is the trivial group {0}.
# The cardinality is 1.
card_H1_B_div_A = 1

# Print the final answer in the format "card(S/A) card(B/S) card(H_1(B/A,Q))"
print(f"{card_S_div_A} {card_B_div_S} {card_H1_B_div_A}")