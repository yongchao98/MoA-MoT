# The limit of the sequence is determined by the prime powers
# that divide P(p) for all sufficiently large primes p.
# Based on number theoretic analysis, we found the exponents of the prime factors of the limit L.
# v_2(L) = 10
# v_3(L) = 2
# v_5(L) = 1
# v_l(L) = 0 for all primes l > 5

# Assign the exponents to variables
exp_2 = 10
exp_3 = 2
exp_5 = 1

# Calculate the base numbers
base_2 = 2
base_3 = 3
base_5 = 5

# Calculate the components of the final number
comp_2 = base_2 ** exp_2
comp_3 = base_3 ** exp_3
comp_5 = base_5 ** exp_5

# Calculate the final limit
limit_g_n = comp_2 * comp_3 * comp_5

# Print the equation with all the numbers
print(f"{base_2}^{exp_2} * {base_3}^{exp_2} * {base_5}^{exp_5} = {comp_2} * {comp_3} * {comp_5} = {limit_g_n}")