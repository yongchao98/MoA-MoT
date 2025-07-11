import numpy as np

def format_polynomial(coeffs):
    """Formats a list of coefficients into a polynomial string."""
    parts = []
    for i, coeff in enumerate(coeffs):
        # Cast to integer if it's a whole number
        if np.isclose(coeff, round(coeff)):
            coeff = int(round(coeff))

        if coeff == 0:
            continue
        
        # Term part
        if i == 0:
            term = str(coeff)
        elif i == 1:
            if coeff == 1:
                term = "x"
            else:
                term = f"{coeff}*x"
        else:
            if coeff == 1:
                term = f"x**{i}"
            else:
                term = f"{coeff}*x**{i}"
        
        parts.append(term)
        
    # Join parts with ' + ' and handle negative signs correctly
    return " + ".join(parts).replace(" + -", " - ")

# Step 3: Define the Poincare polynomials for the component Lie algebras
# P_f4(x) = 1 + 2x + 2x^2 + 2x^3 + x^4 for the filiform algebra f_4
p_f4 = np.array([1, 2, 2, 2, 1])

# P_h3(x) = 1 + 2x + 2x^2 + x^3 for the Heisenberg algebra h_3
p_h3 = np.array([1, 2, 2, 1])

# P_R(x) = 1 + x for the 1D abelian algebra R
p_R = np.array([1, 1])

# Step 4: Compute the final polynomial using the formula
# Numerator = P_f4(x) * P_h3(x)
p_numerator = np.polymul(p_f4, p_h3)

# Final polynomial = Numerator / P_R(x)
p_g, remainder = np.polydiv(p_numerator, p_R)

# The remainder should be zero, indicating exact division
# print(f"Remainder: {remainder}")

# Format and print the final polynomial
final_poly_str = format_polynomial(p_g)
print(final_poly_str)