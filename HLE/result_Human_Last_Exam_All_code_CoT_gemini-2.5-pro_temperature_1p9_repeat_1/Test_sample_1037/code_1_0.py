# The reasoning for the answers is as follows:
# A is the initial object in the category of scales, identified as (id: Z -> Z). So G_A = Z.
# B is the terminal object, which must be (f_0: Z -> {0}) for morphism existence and uniqueness to hold for all scales. So G_B = {0}.
# S is the scale given by the inclusion of Z into the hyperreals, *R. So G_S = *R.

# 1. |S/A| = |G_S / Im(h_AS)| = |*R / Z|
# The cardinality of the hyperreals *R is beth_1. The cardinality of the integers Z is beth_0.
# The cardinality of the quotient |*R / Z| is beth_1 / beth_0 = beth_1.
card_S_div_A = "Beth_1"

# 2. |B/S| = |G_B / Im(h_SB)| = |{0} / Im(h: *R -> {0})|
# The only homomorphism h is the zero map, so its image is {0}.
# The quotient is |{0}/{0}| = 1.
card_B_div_S = 1

# 3. |H_1(B/A, Q)|
# The space B/A is G_B / Im(h_AB) = {0} / Im(h: Z -> {0}), which is {0}/{0}, a single point.
# The first homology group of a point, H_1({pt}, Q), is the trivial group {0}.
# The cardinality of the trivial group is 1.
card_H1_B_div_A = 1

print(f"{card_S_div_A} {card_B_div_S} {card_H1_B_div_A}")