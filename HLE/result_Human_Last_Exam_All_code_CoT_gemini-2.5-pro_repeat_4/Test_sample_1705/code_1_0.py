# Based on the reasoning that the problem is analogous to a simpler piecewise model,
# we found the measure of the set S to be 1/64.

# The number of values in the set is 7.
# In the simplified model g(x) = min(2x, 1), the set of x0 that yields 7 distinct values
# is the interval [1/2^(7-1), 1/2^(7-2)) = [1/64, 1/32).
# The measure (length) of this interval is 1/32 - 1/64 = 1/64.
measure = 1/64

# The problem asks for the measure multiplied by 10^6.
result = measure * 10**6

# The final equation is 1 * 10^6 / 64
print("The final equation is: 1 * 10^6 / 64 = 15625")
print("The Lebesgue measure of S multiplied by 10^6 is:")
print(int(result))