import math

def get_p_adic_valuation(n, p):
    """Computes the p-adic valuation of n."""
    if n == 0:
        return float('inf')
    count = 0
    while n % p == 0:
        n //= p
        count += 1
    return count

# Parameters for the problem
p = 3
k = 3
n = 13

# Apply the formula for the exponent
v_p_n = get_p_adic_valuation(n, p)
floor_term = math.floor(n / (p - 1))
min_term = min(k - 1, floor_term)
exponent = v_p_n + min_term

# The K-group is non-zero if its order is > 1, i.e., exponent > 0
is_nonzero = exponent > 0

print(f"The ring is Z/{p**k}.")
print(f"We are checking the (2n)-th K-group for n = {n}, which is K_{2*n}(\mathbb{{Z}}/{p**k}).")
print(f"The order of this group is given by the formula: p^(v_p(n) + min(k-1, floor(n/(p-1))))")
print(f"For n = {n}, p = {p}, k = {k}:")
print(f"The p-adic valuation v_{p}(n) is v_{p}({n}) = {v_p_n}.")
print(f"The term floor(n/(p-1)) is floor({n}/({p}-1)) = {floor_term}.")
print(f"The term min(k-1, floor(n/(p-1))) is min({k-1}, {floor_term}) = {min_term}.")
print(f"The exponent is {v_p_n} + {min_term} = {exponent}.")
print(f"The order of K_{2*n}(Z/{p**k}) is {p}^{exponent} = {p**exponent}.")
print(f"Since the exponent is greater than 0, the group is non-zero for n = {n}.")
