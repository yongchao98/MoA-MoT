import math

# Set the upper bound for a and b
limit = 1000

# Step 1: Find all primes up to the limit using a sieve
is_prime = [True] * (limit + 1)
is_prime[0] = is_prime[1] = False
for i in range(2, int(math.sqrt(limit)) + 1):
    if is_prime[i]:
        for multiple in range(i * i, limit + 1, i):
            is_prime[multiple] = False
primes = [i for i, p_bool in enumerate(is_prime) if p_bool]

# Step 2: For each prime p, calculate N_p, the number of its powers <= limit.
# A power p^k is counted if k >= 1 and p^k <= limit.
# p^k <= limit  => k * log(p) <= log(limit) => k <= log(limit) / log(p)
n_p_values = []
for p in primes:
    # We only need to consider primes up to the limit itself.
    if p > limit:
        continue
    k_max = int(math.log(limit) / math.log(p))
    if k_max > 0:
        n_p_values.append(k_max)

# Step 3: Calculate T = sum(N_p) and S2 = sum(N_p^2)
# T is the total number of integers > 1 up to the limit that are prime powers.
T = sum(n_p_values)
S2 = sum(x**2 for x in n_p_values)

# Step 4: Calculate the number of un-allowed pairs
# An un-allowed pair is of the form (p^k, q^m) with p != q.
# The number of such pairs is sum_{p!=q} N_p * N_q = (sum N_p)^2 - sum(N_p^2)
num_unallowed_pairs = T**2 - S2

# Step 5: Calculate total pairs and the number of allowed pairs
total_pairs = limit * limit
num_allowed_pairs = total_pairs - num_unallowed_pairs

# Final output: print the equation and the result
print(f"Total pairs (1 <= a,b <= {limit}): {limit} x {limit} = {total_pairs}")
print(f"Number of un-allowed pairs: T^2 - S2 = {T}^2 - {S2} = {num_unallowed_pairs}")
print(f"Number of allowed pairs: {total_pairs} - {num_unallowed_pairs} = {num_allowed_pairs}")
<<<963078>>>