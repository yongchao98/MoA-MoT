import numpy as np
from numpy.polynomial import polynomial as P

def poly_gcd(p1_coeffs, p2_coeffs):
    """
    Computes the greatest common divisor (GCD) of two polynomials
    represented by their coefficient arrays.
    """
    p1 = np.trim_zeros(p1_coeffs, 'b')
    p2 = np.trim_zeros(p2_coeffs, 'b')
    
    while p2.any():
        # Ensure p1 has degree >= p2 for polydiv
        if len(p1) < len(p2):
            p1, p2 = p2, p1
        
        remainder = P.polyrem(p1, p2)
        p1, p2 = p2, np.trim_zeros(remainder, 'b')

    # Normalize the GCD to have a leading coefficient of 1
    if p1.any():
        p1 = p1 / p1[-1]
    
    return p1

# The minimal polynomial of the operator T is p(x) = x^64 + 1
p_coeffs = np.zeros(65)
p_coeffs[64] = 1
p_coeffs[0] = 1

# We check for factors x-1 and x^8-1.
# q(x) = x - 1
q_coeffs = np.array([-1, 1.0])
# r(x) = x^8 - 1
r_coeffs = np.zeros(9)
r_coeffs[8] = 1
r_coeffs[0] = -1

# Calculate the GCDs
gcd_pq = poly_gcd(p_coeffs, q_coeffs)
gcd_pr = poly_gcd(p_coeffs, r_coeffs)

# If the GCD is a non-zero constant (like [1.]), the polynomials are coprime.
# This implies the corresponding linear operators are invertible.
is_coprime1 = (len(gcd_pq) == 1 and gcd_pq[0] != 0)
is_coprime2 = (len(gcd_pr) == 1 and gcd_pr[0] != 0)

dim_ker_T_minus_I = 0 if is_coprime1 else 'non-zero'
dim_ker_T8_minus_I = 0 if is_coprime2 else 'non-zero'

# The dimension of ker(I+...+T^7) is the same as ker(T^8-I)
dim_ker_S8_T = dim_ker_T8_minus_I

print(f"The minimal polynomial of T is p(x) = x^64 + 1.")
print(f"The GCD of p(x) and x-1 is {gcd_pq}. This means dim(ker(T-I)) = {dim_ker_T_minus_I}.")
print(f"The GCD of p(x) and x^8-1 is {gcd_pr}. This means dim(ker(T^8-I)) = {dim_ker_T8_minus_I}.")
print()

# Final calculation
# dim H^2(G, M) = dim H^1(G, M) = dim(ker(T-I)) + dim(ker(I+...+T^7))
result = dim_ker_T_minus_I + dim_ker_S8_T
print("The dimension of the cohomology group is calculated as:")
print(f"dim H^2(G, M) = dim(ker(T-I)) + dim(ker(I+T+...+T^7))")
print(f"dim H^2(G, M) = {dim_ker_T_minus_I} + {dim_ker_S8_T} = {result}")