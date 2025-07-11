import numpy as np

def print_polynomial(coeffs):
    """
    Prints a polynomial from a list of coefficients.
    Example: [1, 2, 3] -> 1 + 2*x + 3*x**2
    """
    poly_str = []
    for i, c in enumerate(coeffs):
        if c == 0:
            continue
        
        # Handle coefficient part
        if c == 1 and i > 0:
            coeff_part = ""
        else:
            coeff_part = str(c)

        if i > 0:
            # Handle variable part
            if i == 1:
                var_part = "x"
            else:
                var_part = f"x**{i}"
            
            # Add * if coefficient is present
            if coeff_part:
                 term = f"{coeff_part}*{var_part}"
            else:
                 term = var_part
        else:
            term = coeff_part

        poly_str.append(term)
    
    # Join with ' + ' and replace '+ -' with '- ' for negative coefficients
    final_str = " + ".join(poly_str).replace("+ -", "- ")
    print(final_str)

# The Poincaré polynomial of g is given by the product of three polynomials:
# P_g(x) = (1+x) * (1+x+x^2+x^3) * (1+x+x^2)

p1_coeffs = [1, 1]  # Represents 1+x
p2_coeffs = [1, 1, 1, 1]  # Represents 1+x+x^2+x^3
p3_coeffs = [1, 1, 1]  # Represents 1+x+x^2

# Multiply the polynomials
# Result1 = p1 * p2
result1_coeffs = np.polymul(p1_coeffs, p2_coeffs)

# Result = Result1 * p3. We reverse lists for numpy convention (highest degree first).
final_coeffs_np = np.polymul(result1_coeffs[::-1], p3_coeffs[::-1])

# Reverse back to lowest degree first for printing
final_coeffs = final_coeffs_np.tolist()[::-1]

# Print the final resulting polynomial equation
print("The Poincaré polynomial is P(x) =")
print_polynomial(final_coeffs)