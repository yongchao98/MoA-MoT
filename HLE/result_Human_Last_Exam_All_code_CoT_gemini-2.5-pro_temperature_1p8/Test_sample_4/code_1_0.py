import numpy as np

def poly_to_string(coeffs, var='x'):
    """Converts a list of polynomial coefficients to a string."""
    parts = []
    for i, c in enumerate(coeffs):
        if c == 0:
            continue
        # Term part
        if i == 0:
            term = str(c)
        elif i == 1:
            term = f"{c}{var}" if c != 1 else var
        else:
            term = f"{c}{var}^{i}" if c != 1 else f"{var}^{i}"
        
        # Sign part
        if not parts:
             parts.append(term)
        elif c > 0:
            parts.append(f" + {term}")
        else:
            # this case won't happen with the current problem
            parts.append(f" - {str(abs(c)) + term if abs(c)!=1 else term}")

    return "".join(parts[::-1]) # reverse for standard poly order

# Coefficients of the polynomials, from degree 0 to max
# P_n4(x) = 1 + 2x + 2x^2 + 2x^3 + x^4
p_n4 = [1, 2, 2, 2, 1]
# P_h3(x) / P_R(x) = (1 + 2x + 2x^2 + x^3) / (1+x) = 1 + x + x^2
p_quotient = [1, 1, 1]

# Perform polynomial multiplication
p_g = np.polymul(p_n4[::-1], p_quotient[::-1])[::-1]
p_g_coeffs = [int(c) for c in p_g]

final_poly_str = poly_to_string(p_g_coeffs, 'x')

print("The Poincaré polynomial P(x) is the product of P_n4(x) and (P_h3(x)/(1+x)).")
print(f"P_n4(x) = {poly_to_string([1,2,2,2,1])}")
print(f"P_h3(x) / (1+x) = {poly_to_string([1,1,1])}")
print(f"\nP(x) = ({poly_to_string([1,2,2,2,1])}) * ({poly_to_string([1,1,1])})")
print(f"P(x) = {final_poly_str}")

# Final Answer formatting as requested
equation_parts = []
for i, c in enumerate(p_g_coeffs):
    if c == 0:
        continue
    
    if i == 0:
        part_str = f"{c}"
    elif i == 1:
        part_str = f"{c}x" if c != 1 else "x"
    else:
        part_str = f"{c}x^{i}" if c != 1 else f"x^{i}"
    equation_parts.append(part_str)

final_equation = " + ".join(equation_parts[::-1])
print(f"\nThe final equation is:")
print(f"{final_equation}")

# Final Answer Wrapper
print("\nFinal Answer in specified format:")
print(f"<<<{final_equation}>>>")