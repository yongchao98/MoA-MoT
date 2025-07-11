import numpy as np

# Coefficients of Q1(x) = 1 + x + x^2 + x^3
# for a Jordan block of size 3
Q1_coeffs = [1, 1, 1, 1]

# Coefficients of Q2(x) = 1 + x + x^2
# for a Jordan block of size 2
Q2_coeffs = [1, 1, 1]

# Compute coefficients of Q(x) = Q1(x) * Q2(x)
Q_coeffs = np.polymul(Q1_coeffs, Q2_coeffs)

# Compute coefficients of P(x) = (1+x) * Q(x)
one_plus_x = [1, 1]
P_coeffs = np.polymul(one_plus_x, Q_coeffs)

# Format the output string for the polynomial
poly_str = []
for i, coeff in enumerate(P_coeffs):
    coeff = int(coeff)
    if coeff == 0:
        continue
    # Term part
    if i == 0:
        term_str = str(coeff)
    elif i == 1:
        if coeff == 1:
            term_str = "x"
        else:
            term_str = f"{coeff}x"
    else:
        if coeff == 1:
             term_str = f"x^{i}"
        else:
            term_str = f"{coeff}x^{i}"
    poly_str.append(term_str)

print("The Poincare polynomial is P(x) = " + " + ".join(poly_str))