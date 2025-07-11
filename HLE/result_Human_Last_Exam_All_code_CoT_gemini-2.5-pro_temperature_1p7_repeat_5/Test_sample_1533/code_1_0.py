# The side lengths of the triangle ABC are denoted by a, b, and c.
# a is the length of side BC, opposite to vertex A.
# b is the length of side AC, opposite to vertex B.
# c is the length of side AB, opposite to vertex C.

# Based on the geometric derivation using the Incenter-Excenter Lemma and Ptolemy's Theorem,
# the ratio of the lengths of the segments BM and MI is found to be (a + c) / b.

# To fulfill the requirement of showing each number in the final equation,
# we can write the expression with explicit coefficients.
# The equation is BM/MI = (1*a + 1*c) / (1*b).

coef_a = 1
coef_c = 1
coef_b = 1

# Print the final result as a formatted string.
print("The ratio BM/MI is expressed in terms of the side lengths a, b, and c.")
print("The final equation is:")
print(f"BM / MI = ({coef_a}*a + {coef_c}*c) / {coef_b}*b")
