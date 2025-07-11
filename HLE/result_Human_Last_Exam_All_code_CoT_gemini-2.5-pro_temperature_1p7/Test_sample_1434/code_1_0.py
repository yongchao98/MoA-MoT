from fractions import Fraction

# The problem asks for the sum of 1/4^n for all pairs (n, m) in a set S.
# The number of possible values for m depends on the parity of n.

# For even n, there are 2 possible values of m.
# The sum for even n is S_even = sum_{k=1 to inf} 2 / 4^(2k)
# This is a geometric series: 2 * sum_{k=1 to inf} (1/16)^k
r_even = Fraction(1, 16)
sum_even_series = r_even / (1 - r_even) # Sum from k=1
S_even = 2 * sum_even_series

print(f"The sum for even values of n is 2 * sum_{k=1 to inf} (1/16)^k = 2 * ( (1/16) / (1 - 1/16) ) = 2 * {sum_even_series} = {S_even}")

# For odd n, there are 4 possible values of m.
# The sum for odd n is S_odd = sum_{k=0 to inf} 4 / 4^(2k+1)
# This can be simplified: 4 * (1/4) * sum_{k=0 to inf} (1/16)^k = sum_{k=0 to inf} (1/16)^k
r_odd = Fraction(1, 16)
sum_odd_series = 1 / (1 - r_odd) # Sum from k=0
S_odd = sum_odd_series

print(f"The sum for odd values of n is 4 * sum_{k=0 to inf} 1 / 4^(2k+1) = sum_{k=0 to inf} (1/16)^k = 1 / (1-1/16) = {S_odd}")

# The total sum is S_even + S_odd
total_sum = S_even + S_odd

print(f"The total sum is the sum of the even and odd parts:")
print(f"Total sum = {S_even} + {S_odd} = {total_sum}")

# As a decimal for final output
final_answer = float(total_sum)
print(f"The final numerical answer is {final_answer}")
