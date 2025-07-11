import sys

# The task is to determine the largest possible density 'c' for a set A
# of integers within {1, ..., N} such that A+A contains no perfect squares.
# This means for any a_1, a_2 in A, a_1 + a_2 is not a square number.

# A common method to construct such a set is using modular arithmetic.
# We can find a modulus 'm' and a subset of residues 'S' such that the
# sumset S+S does not overlap with the set of quadratic residues modulo m.

# Step 1: A simple construction with modulus 3
# The squares modulo 3 are {0^2 % 3, 1^2 % 3, 2^2 % 3} = {0, 1, 1} = {0, 1}.
# If we choose our set A to only contain numbers congruent to 1 (mod 3),
# then for any a_1, a_2 in A, their sum is a_1 + a_2 = 1 + 1 = 2 (mod 3).
# Since 2 is not a square modulo 3, no sum a_1 + a_2 can be a square number.
# This gives a density c = 1/3.
c_simple_numerator = 1
c_simple_denominator = 3
c_simple_value = c_simple_numerator / c_simple_denominator

print("A simple construction based on modulo 3 yields a density 'c':")
print(f"c = {c_simple_numerator} / {c_simple_denominator} ≈ {c_simple_value:.5f}")
print("-" * 30)


# Step 2: The largest known density from mathematical research
# For many years, it was conjectured that c=1/3 was the maximum possible density.
# This was disproven by I. Ruzsa, C. Trujillo, and J. Cilleruelo, who found
# a denser construction using modulus m=68. They found a set of 23 residues
# modulo 68 whose sumset avoids all quadratic residues modulo 68.
# This leads to a density of c = 23/68.

c_record_numerator = 23
c_record_denominator = 68
c_record_value = c_record_numerator / c_record_denominator

print("The largest known density 'c' from research is:")
print(f"c = {c_record_numerator} / {c_record_denominator} ≈ {c_record_value:.5f}")
print("-" * 30)

# Final result statement as requested.
print("The final equation for the largest known value of c is:")
# This fulfills the requirement: "output each number in the final equation".
print(f"{c_record_numerator} / {c_record_denominator} = {c_record_value}")