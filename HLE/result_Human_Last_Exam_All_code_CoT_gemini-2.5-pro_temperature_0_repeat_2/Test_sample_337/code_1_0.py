# Number of prisoners
n = 15

# The total number of possible hat configurations is 2^n, since each of the n
# prisoners can have one of two hat colors (black or white).
total_configurations = 2**n

# The optimal strategy involves the prisoners agreeing on a "target parity" for the
# number of black hats (e.g., they will assume the total is even). This strategy
# allows them to win for exactly half of all possible configurations.
# The number of winning configurations is 2^(n-1).
winning_configurations = 2**(n - 1)

# The maximal probability is the ratio of the number of winning configurations
# to the total number of configurations.
probability = winning_configurations / total_configurations

# Print the explanation and the final calculation, showing each number.
print(f"The number of prisoners is {n}.")
print(f"The total number of possible hat configurations is 2^{n}, which is {total_configurations}.")
print(f"An optimal strategy allows the prisoners to win on half of these configurations.")
print(f"The number of winning configurations is 2^({n}-1), which is {winning_configurations}.")
print(f"The maximal probability of being released is the ratio of winning configurations to the total configurations:")
print(f"P_max = {winning_configurations} / {total_configurations} = {probability}")