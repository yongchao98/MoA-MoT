# Let a, b, and c be the side lengths of triangle ABC,
# where a is the length of side BC, b is the length of side CA,
# and c is the length of side AB.

# We define these as string variables to construct the final formula.
side_a = 'a'
side_b = 'b'
side_c = 'c'

# From the geometric derivation, the ratio BM/MI is found to be (a+c)/b.
# We will print this formula. The components of the equation are:
# Numerator part 1: a
# Numerator part 2: c
# Denominator: b

print("The ratio BM / MI in terms of the side lengths a, b, and c is:")
print(f"BM / MI = ({side_a} + {side_c}) / {side_b}")
