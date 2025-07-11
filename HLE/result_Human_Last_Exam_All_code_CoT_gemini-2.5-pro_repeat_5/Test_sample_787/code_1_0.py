# The exponents of the prime factors of the limit L are determined by the
# analysis above.
# a_2 = v_2(L) = 10
# a_3 = v_3(L) = 2
# a_5 = v_5(L) = 1
# a_q = v_q(L) = 0 for all other primes q.

# Define the exponents
a2 = 10
a3 = 2
a5 = 1

# Calculate the limit L
limit_L = (2**a2) * (3**a3) * (5**a5)

# Print the final equation and the result
print(f"The limit is L = 2^{a2} * 3^{a3} * 5^{a5}")
print(f"L = {2**a2} * {3**a3} * {5**a5}")
print(f"L = {limit_L}")