# The problem asks for the smallest preference profile size n for which
# a committee of size k=100 can exist that satisfies a given representation
# property (PJR/EJR) while leaving a specific voter (voter 1) unsatisfied.

# Let's define the committee size.
k = 100

# --- Part 1: Proportional Justified Representation (s1) ---

# We need a committee W that satisfies PJR but for which voter 1 is unsatisfied.
# A voter 'i' being unsatisfied means their approval ballot A(i) has no candidates
# in the committee W, i.e., |A(i) intersect W| = 0.

# The PJR property is violated if there exists a group of voters N' of size l,
# where l >= n/k (n is the total number of voters), that are cohesive (all
# approve at least one common candidate), and are all unsatisfied.

# Consider a group N' consisting only of voter 1. The size is l=1.
# If n <= k, then n/k <= 1, so l=1 >= n/k.
# This means voter 1 forms a valid group for the PJR check.
# If we have a committee W where voter 1 is unsatisfied, this group N'={voter 1}
# is cohesive and all its members (just voter 1) are unsatisfied. This is a
# PJR violation.
# Therefore, if n <= k, it's impossible for a committee to satisfy PJR and
# leave voter 1 unsatisfied.

# To make it possible, we must ensure the group {voter 1} is not large enough
# to trigger a PJR check. We need its size l=1 to be smaller than n/k.
# This requires: 1 < n/k, which is equivalent to n > k.
# The smallest integer n satisfying n > k is k + 1.
# For such a profile to exist, n = k+1 is the minimum size.

s1_k_val = k
s1_1_val = 1
s1 = s1_k_val + s1_1_val

# --- Part 2: Extended Justified Representation (s2) ---

# EJR is a stronger property. An EJR violation occurs if there is a group N'
# (size l >= n/k) and an integer j >= 1 such that all voters in N' unanimously
# approve a set of candidates C with |C| >= j, but for every voter v in N',
# their representation is less than j (i.e., |A(v) intersect W| < j).

# We need a committee W that satisfies EJR but leaves voter 1 unsatisfied.
# Voter 1 is unsatisfied, so |A(1) intersect W| = 0. This implies that no
# candidate from voter 1's ballot, {a, b, c, x}, can be in W.

# Consider the group N' = {voter 1, voter 2, ..., voter 6}. Its size is m = 6.
# Let's analyze this group:
# - All 6 voters approve the set C = {a, b, c}. So |C| = 3.
# - This group is therefore "j-cohesive" for j=1, 2, and 3.

# Let's check the representation of these voters for j=2.
# Since {a, b, c} are not in W, the representation of each voter in this group is at most 1
# (e.g., voter 2 is 1-represented if y is in W, but not more).
# So, for every voter v in this group N', |A(v) intersect W| <= 1, which is < 2.

# This group N' will cause an EJR violation for j=2 IF the group is large enough.
# A violation occurs if |N'| >= n/k.
# |N'| = m = 6. The violation happens if 6 >= n/k, or n <= 6k.
# To prevent this violation and allow for a committee that satisfies EJR,
# the group must not be large enough. We need |N'| < n/k.
# This requires: m < n/k  --->  6 < n/100  --->  600 < n.
# The smallest integer n greater than 600 is 601.

m_val = 6
s2_m_val = m_val
s2_k_val = k
s2_1_val = 1
s2 = s2_m_val * s2_k_val + s2_1_val

# Print the calculations and the final answer.
print(f"The committee size is k = {k}")
print(f"The size of the core group for the EJR check is m = {m_val}")
print()
print("Calculation for s1 (PJR):")
print(f"s1 = k + 1 = {s1_k_val} + {s1_1_val} = {s1}")
print()
print("Calculation for s2 (EJR):")
print(f"s2 = m * k + 1 = {s2_m_val} * {s2_k_val} + {s2_1_val} = {s2}")
print()
print(f"The final solution as a pair (s1, s2) is: ({s1}, {s2})")