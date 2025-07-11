# The limit of the sequence g_n is the greatest common divisor (GCD)
# of P(p) for all sufficiently large prime numbers p. Let this limit be L.
# L = lim_{n->inf} gcd_{p > n, p is prime} P(p).
# Based on mathematical analysis of the divisibility of the polynomial P(X),
# the prime factorization of L can be determined by finding the minimum
# q-adic valuation of P(p) over all large primes p.

# The analysis yields the following prime powers for L:
# 1. Power of 2: The minimum 2-adic valuation v_2(P(p)) for p > 2 is 10.
# 2. Power of 3: The minimum 3-adic valuation v_3(P(p)) for p > 3 is 2.
# 3. Power of 5: The minimum 5-adic valuation v_5(P(p)) for p > 5 is 1.
# 4. Powers of other primes: For any prime q > 5, the minimum valuation is 0.

# Therefore, L = 2^10 * 3^2 * 5^1. We calculate this value below.

# Exponents for the prime factorization of the limit
v2 = 10
v3 = 2
v5 = 1

# Bases of the prime factors
base2 = 2
base3 = 3
base5 = 5

# Calculate the limit by multiplying the prime powers
limit = (base2 ** v2) * (base3 ** v3) * (base5 ** v5)

# Output the result in the specified format
print(f"{base2}^{v2} * {base3}^{v3} * {base5}^{v5} = {limit}")