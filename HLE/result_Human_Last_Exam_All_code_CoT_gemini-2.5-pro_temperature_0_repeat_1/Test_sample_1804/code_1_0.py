def poly_div_rem(N, D, p):
    """
    Performs polynomial long division over a finite field F_p.
    N and D are lists of coefficients in descending order of power.
    Returns the remainder as a list of coefficients.
    """
    num = list(N)
    den = list(D)
    
    if len(den) > len(num):
        return num

    # Modular inverse of the leading coefficient of the divisor
    try:
        inv = pow(den[0], p - 2, p)
    except (ValueError, ZeroDivisionError):
        # This case should not be reached with monic divisors
        return None

    # Perform division steps
    while len(num) >= len(den):
        q_coeff = (num[0] * inv) % p
        
        # Subtract the multiple of the divisor from the dividend
        for i in range(len(den)):
            sub_term = (q_coeff * den[i]) % p
            num[i] = (num[i] - sub_term + p) % p
        
        # Remove the leading zero
        num.pop(0)

    # Clean up any remaining leading zeros in the remainder
    while num and num[0] == 0:
        num.pop(0)
        
    return num

def is_irreducible(a, p=7):
    """
    Checks if the polynomial x^5 + ax + 3 is irreducible over the finite field F_p.
    """
    poly_coeffs = [1, 0, 0, 0, a, 3]

    # 1. Check for roots in F_p. If a root exists, the polynomial is reducible.
    for c in range(p):
        val = (pow(c, 5, p) + (a * c) % p + 3) % p
        if val == 0:
            return False

    # 2. Check for irreducible quadratic factors.
    # A degree 5 polynomial is reducible if it has a factor of degree 1 or 2.
    # We already checked for degree 1 (roots). Now check for degree 2.
    
    # Iterate through all monic quadratics x^2 + bx + c
    for b in range(p):
        for c in range(p):
            # First, check if the quadratic x^2 + bx + c is itself irreducible.
            # It is irreducible if it has no roots in F_p.
            is_quad_reducible = False
            for x_val in range(p):
                if (x_val*x_val + b*x_val + c) % p == 0:
                    is_quad_reducible = True
                    break
            
            # If the quadratic is reducible, we skip it.
            if is_quad_reducible:
                continue

            # If the quadratic is irreducible, check if it divides our polynomial.
            divisor_coeffs = [1, b, c]
            remainder = poly_div_rem(poly_coeffs, divisor_coeffs, p)
            
            # If the remainder is zero, our polynomial is divisible and thus reducible.
            if not remainder:
                return False

    # If no roots and no irreducible quadratic factors are found, it's irreducible.
    return True

# Main part of the script
F_p = 7
A = []
for a_val in range(F_p):
    if is_irreducible(a_val, F_p):
        A.append(a_val)

# Calculate the required values from the set A
min_A = min(A)
max_A = max(A)
len_A = len(A)

# Perform the final calculation
result = max_A**min_A - len_A

# Print the components of the final equation as requested
print(f"The set A is: {A}")
print(f"The value for max(A) in the equation is: {max_A}")
print(f"The value for min(A) in the equation is: {min_A}")
print(f"The value for |A| in the equation is: {len_A}")
print(f"The final result of {max_A}^{min_A} - {len_A} is: {result}")