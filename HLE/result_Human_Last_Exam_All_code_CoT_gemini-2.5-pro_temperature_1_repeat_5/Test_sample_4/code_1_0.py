# The Betti numbers for the Lie algebra g.
betti_numbers = [1, 3, 6, 8, 6, 3, 1]

# Constructing the Poincaré polynomial string representation.
polynomial_terms = []
for i, b in enumerate(betti_numbers):
    if b == 0:
        continue
    if i == 0:
        term = str(b)
    elif i == 1:
        if b == 1:
            term = "x"
        else:
            term = f"{b}*x"
    else:
        if b == 1:
            term = f"x^{i}"
        else:
            term = f"{b}*x^{i}"
    polynomial_terms.append(term)

polynomial_string = " + ".join(polynomial_terms)

print("The Poincaré polynomial is:")
print(polynomial_string)

# The final answer asked to output each number in the final equation.
# Let's format it cleanly.
equation = "P(x) = "
for i, b in enumerate(betti_numbers):
    if b > 0:
        if i == 0:
            equation += f"{b}"
        elif i == 1:
            equation += f" + {b}*x"
        else:
            equation += f" + {b}*x^{i}"
print("\nFormatted equation:")
print(equation)

# Just to be sure, let's output the numbers directly.
print("\nThe coefficients (Betti numbers b_0 to b_6) are:")
for number in betti_numbers:
    print(number)
