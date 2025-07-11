import math

# An internal adjunction from [m] to [n] in the simplex category Delta
# corresponds to a monotone map L: [m] -> [n] where L(0) = 0.
# Here, m=23 and n=37.

m = 23
n = 37

# The condition L(0)=0 means we need to count the number of non-decreasing sequences
# L(1), L(2), ..., L(23) with values chosen from {0, 1, ..., 37}.

# This is a combination with repetition problem (stars and bars).
# The number of values to choose is k = m.
k = m

# The number of options for each value is the size of the set {0, 1, ..., n}.
num_options = n + 1

# The formula for combinations with repetition is C(n' + k - 1, k),
# where n' is the number of options.
n_comb = num_options + k - 1
k_comb = k

# Calculate the binomial coefficient C(60, 23)
result = math.comb(n_comb, k_comb)

print(f"The number of internal adjunctions from [{m}] to [{n}] is given by the binomial coefficient C(n', k), where:")
print(f"n' is the number of items to choose from with repetition, and k is the number of items being chosen.")
print(f"k = m = {m}")
print(f"n' = (n + 1) = {n} + 1 = {num_options}")
print(f"The number of adjunctions is C(n' + k - 1, k) = C({num_options} + {k} - 1, {k}) = C({n_comb}, {k_comb}).")
print(f"C({n_comb}, {k_comb}) = {result}")