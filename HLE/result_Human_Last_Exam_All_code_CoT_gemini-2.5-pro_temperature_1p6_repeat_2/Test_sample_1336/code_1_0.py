# The simple group in question is S = PSL(2, p) for a prime p > 5.
# The number of smooth coverings is equal to the number of non-isomorphic covering groups of S.
# This number is, in turn, equal to the number of subgroups of the Schur multiplier M(S).

# For S = PSL(2, p) with p > 5 prime, the Schur multiplier M(S) is the cyclic group of order 2, C_2.
# We need to count the number of subgroups of C_2.

# A cyclic group of order n has one subgroup for each divisor of n.
# For C_2, n=2. The divisors of 2 are 1 and 2. So there are two subgroups.

# Subgroup 1: The trivial subgroup {e}. This corresponds to one covering group (PSL(2,p)).
covering_from_trivial_subgroup = 1

# Subgroup 2: The group C_2 itself. This corresponds to another covering group (SL(2,p)).
covering_from_C2_subgroup = 1

# The total number of coverings is the sum.
total_coverings = covering_from_trivial_subgroup + covering_from_C2_subgroup

# The problem asks to output the numbers in the final equation.
print(f"The number of coverings from the two subgroups of the Schur multiplier are {covering_from_trivial_subgroup} and {covering_from_C2_subgroup}.")
print(f"The total number of smooth coverings is the sum: {covering_from_trivial_subgroup} + {covering_from_C2_subgroup} = {total_coverings}")