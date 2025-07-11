# The side length of the square table.
N = 2024

# Based on the analysis, the minimum number of diamonds required to form a "frozen"
# arrangement, where no diamond can be legally moved to an adjacent cell, is N.
# An example of such a minimal frozen arrangement is placing diamonds along the main diagonal.
m_min = N

# The question asks for the largest number of diamonds, k, such that ANY arrangement
# of k diamonds is NOT frozen. This means k must be one less than the minimum
# size of a frozen arrangement.
k = m_min - 1

# Printing the final equation and the result.
print(f"The side of the square is N = {N}.")
print(f"The minimum number of diamonds in a frozen arrangement is m_min = {m_min}.")
print("The largest value k such that any arrangement of k diamonds is not frozen is k = m_min - 1.")
print(f"The final equation is: {m_min} - 1 = {k}")