import numpy as np

# P_ker_A(x) = 1 + x + x^2 + x^3
p_ker_A = np.array([1, 1, 1, 1])

# P_ker_B(x) = 1 + x + x^2
p_ker_B = np.array([1, 1, 1])

# P_ker(x) = P_ker_A(x) * P_ker_B(x)
p_ker = np.polymul(p_ker_A, p_ker_B)

# P_g(x) = (1 + x) * P_ker(x)
p_g = np.polymul(np.array([1, 1]), p_ker)

# Format the output string for the polynomial
# The coefficients are ordered from highest degree to lowest
coeffs = p_g[::-1]
poly_str = []
for i, c in enumerate(coeffs):
    if c == 0:
        continue
    # Do not print coefficient if it is 1 and not the constant term
    coeff_str = "" if (c == 1 and i > 0) else str(int(c))
    
    # Do not print power if it is 0
    if i == 0:
        var_str = ""
    # Do not print power if it is 1
    elif i == 1:
        var_str = "x"
    else:
        var_str = f"x^{i}"
    
    # Combine coefficient and variable part
    term = f"{coeff_str}{var_str}"
    
    poly_str.append(term)

print("The Poincaré polynomial is: P(x) =", " + ".join(poly_str))
# The final task asks to output each number in the final equation.
# This means we should print the polynomial in its expanded form.
# The code above already generates the string representation.