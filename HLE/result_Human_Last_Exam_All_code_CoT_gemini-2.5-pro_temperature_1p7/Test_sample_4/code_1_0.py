# The Betti numbers for the Lie algebra g.
b = [1, 3, 6, 8, 6, 3, 1]

# Construct the polynomial string
poly_str = []
for i, coeff in enumerate(b):
    if coeff == 0:
        continue
    if i == 0:
        poly_str.append(f"{coeff}")
    elif i == 1:
        if coeff == 1:
            poly_str.append("x")
        else:
            poly_str.append(f"{coeff}*x")
    else:
        if coeff == 1:
            poly_str.append(f"x^{i}")
        else:
            poly_str.append(f"{coeff}*x^{i}")

final_poly = " + ".join(poly_str)

# The question asks to output the equation with each number.
# Let's create a formatted string for the final answer.
# P(x) = b_0 + b_1*x + b_2*x^2 + ...
equation_parts = []
equation_parts.append(str(b[0]))
equation_parts.append(f"{b[1]}*x")
for i in range(2, len(b)):
    equation_parts.append(f"{b[i]}*x^{i}")

final_equation_str = " + ".join(equation_parts)

print("The Poincaré polynomial is:")
print(f"P(x) = {final_equation_str}")