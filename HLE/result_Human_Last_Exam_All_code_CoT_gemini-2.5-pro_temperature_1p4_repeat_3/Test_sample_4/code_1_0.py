import numpy as np

def poly_mult(p1_coeffs, p2_coeffs):
    """Multiplies two polynomials given as lists of coefficients."""
    return np.convolve(p1_coeffs, p2_coeffs).tolist()

def format_poly(p_coeffs, name="P(x)"):
    """Formats a polynomial for printing."""
    terms = []
    # Iterate in reverse to get descending powers of x
    for i, coeff in enumerate(p_coeffs):
        coeff = int(round(coeff))
        if coeff == 0:
            continue
        
        # Format coefficient
        if coeff == 1 and i > 0:
            c_str = ""
        else:
            c_str = str(coeff)

        # Format variable part
        if i == 0:
            v_str = ""
        elif i == 1:
            v_str = "x"
        else:
            v_str = f"x^{i}"
        
        if c_str and v_str:
            terms.append(f"{c_str}{v_str}")
        else:
            terms.append(f"{c_str or v_str}")

    if not terms:
        return f"{name} = 0"
    return f"{name} = {' + '.join(terms)}"

print("Step 1: Define the factors of the Poincaré polynomial for g.")

p1 = {'name': '1+x', 'coeffs': [1, 1]}
p2 = {'name': '1+x+x^2', 'coeffs': [1, 1, 1]}
p3 = {'name': '1+x+x^2+x^3', 'coeffs': [1, 1, 1, 1]}

print("The Poincaré polynomial is a product of three factors:")
print(format_poly(p1['coeffs'], p1['name']))
print(format_poly(p2['coeffs'], p2['name']))
print(format_poly(p3['coeffs'], p3['name']))

print("\nStep 2: Compute the product of these polynomials.")

print("First, multiply (1+x+x^2) and (1+x+x^2+x^3).")
prod1_coeffs = poly_mult(p2['coeffs'], p3['coeffs'])
print(format_poly(prod1_coeffs, f"({p2['name']})({p3['name']})"))

print("\nNow, multiply the result by (1+x).")
final_coeffs = poly_mult(prod1_coeffs, p1['coeffs'])

print("\nThe final Poincaré polynomial for g is:")

# Custom formatting for the final answer
terms = []
for i, coeff in enumerate(final_coeffs):
    coeff = int(round(coeff))
    if coeff == 0:
        continue
    c_str = str(coeff)
    if i == 0:
        v_str = ""
    elif i == 1:
        v_str = "x"
    else:
        v_str = f"x^{i}"
    
    if v_str:
        if coeff == 1:
            terms.append(v_str)
        else:
            terms.append(f"{c_str} {v_str}")
    else:
        terms.append(c_str)

final_poly_str = " + ".join(terms)
print(f"P_g(x) = {final_poly_str}")

# Storing the answer in the requested format
answer = f"1 + 3x + 5x^2 + 6x^3 + 5x^4 + 3x^5 + x^6"