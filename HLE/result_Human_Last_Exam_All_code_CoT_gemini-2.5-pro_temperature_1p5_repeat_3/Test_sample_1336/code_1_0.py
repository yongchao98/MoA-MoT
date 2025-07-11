# The number of smooth coverings corresponds to the number of quasi-simple covering
# groups of S = PSL(2, p). This is determined by the Schur Multiplier of S.

# For p > 5 prime, the Schur Multiplier M(PSL(2, p)) is the cyclic group of order 2.
order_of_schur_multiplier = 2

# The number of covering groups is equal to the number of subgroups of the Schur Multiplier.
# For a cyclic group of order n, this is the number of divisors of n, tau(n).
# We calculate tau(2).

n = order_of_schur_multiplier
num_divisors = 0
divisors = []
for i in range(1, n + 1):
    if n % i == 0:
        num_divisors += 1
        divisors.append(i)

# The total number of smooth coverings is the result of this calculation.
total_coverings = num_divisors

# As requested, we output the final equation.
# The equation shows that the number of coverings is tau(2).
print(f"The number of coverings = tau(order of M(PSL(2,p)))")
print(f"tau({order_of_schur_multiplier}) = {total_coverings} (since the divisors of {order_of_schur_multiplier} are {divisors})")
print(f"Thus, the total number of such smooth coverings is {total_coverings}.")
