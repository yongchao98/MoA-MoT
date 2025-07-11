#
# This script calculates the number of multiplication operations in the fully expanded
# boolean expression for the s2 sum bit of a 3-bit binary adder.
#
# A binary number A is represented by a2, a1, a0.
# A binary number B is represented by b2, b1, b0.
# The sum S is represented by s3, s2, s1, s0.
#
# The formula for s2 is:
# s2 = (a2' * b2' * c2) + (a2' * b2 * c2') + (a2 * b2' * c2') + (a2 * b2 * c2)
#

# --- Analysis of the first component: (a2' * b2' * c2) ---
# The expression for c2 is (a1*b1) + (a1*a0*b0) + (b1*a0*b0).
# The terms of c2 have 2, 3, and 3 literals, respectively.
# We multiply each by (a2'*b2'), which has 2 literals.
# Resulting terms' literals: (2+2), (3+2), (3+2) => 4, 5, 5
# Multiplications for a term with N literals is N-1.
muls_part1 = (4 - 1) + (5 - 1) + (5 - 1)
print(f"Number of multiplications in (a2' * b2' * c2) expanded part: {muls_part1}")

# --- Analysis of the second component: (a2' * b2 * c2') ---
# The expression for c2' has 5 terms, each with 2 literals.
# We multiply each by (a2'*b2), which has 2 literals.
# Resulting terms' literals: 5 terms, each with (2+2) = 4 literals.
# Multiplications: 5 * (4-1)
muls_part2 = 5 * (4 - 1)
print(f"Number of multiplications in (a2' * b2 * c2') expanded part: {muls_part2}")

# --- Analysis of the third component: (a2 * b2' * c2') ---
# This is symmetric to the second component.
# We multiply each of the 5 terms of c2' (2 literals each) by (a2*b2') (2 literals).
# Resulting terms' literals: 5 terms, each with (2+2) = 4 literals.
# Multiplications: 5 * (4-1)
muls_part3 = 5 * (4 - 1)
print(f"Number of multiplications in (a2 * b2' * c2') expanded part: {muls_part3}")

# --- Analysis of the fourth component: (a2 * b2 * c2) ---
# This is symmetric to the first component.
# We multiply each term of c2 by (a2*b2).
# Resulting terms' literals: (2+2), (3+2), (3+2) => 4, 5, 5
# Multiplications: (4-1) + (5-1) + (5-1)
muls_part4 = (4 - 1) + (5 - 1) + (5 - 1)
print(f"Number of multiplications in (a2 * b2 * c2) expanded part: {muls_part4}")

# --- Total Calculation ---
total_multiplications = muls_part1 + muls_part2 + muls_part3 + muls_part4
print(f"\nThe total number of multiplication operations is the sum of the parts:")
print(f"{muls_part1} + {muls_part2} + {muls_part3} + {muls_part4} = {total_multiplications}")

print("\n<<<52>>>")