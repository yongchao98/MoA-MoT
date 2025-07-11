import sympy

# Define 'k' as a symbol for our polynomial variable
k = sympy.Symbol('k')

# Step 1: Define the chromatic polynomials for K3 and K4
P_K3 = k * (k - 1) * (k - 2)
P_K4 = k * (k - 1) * (k - 2) * (k - 3)

# Step 2: Calculate the chromatic polynomial of the given graph G
# This is based on the decomposition explained above:
# P(G, k) = P(K3, k) * (k-1) + P(K4, k) * (k-2)
P_G = P_K3 * (k - 1) + P_K4 * (k - 2)

# Step 3: Expand the polynomial to get the standard form
expanded_P_G = sympy.expand(P_G)

# Step 4: Format the output to clearly show the equation
# with all its numerical coefficients.
# We will represent the polynomial as a string.
poly_str = str(expanded_P_G).replace('**', '^').replace('*', ' * ')

# Print the final result in a readable equation format
print("The chromatic polynomial of the graph is P(k):")
print(f"P(k) = {poly_str}")

# We can also iterate through the coefficients and print them term by term
print("\nIn expanded form, the equation with each number is:")
terms = []
# Get coefficients from the expanded polynomial, in descending order of power
coeffs = sympy.Poly(expanded_P_G, k).all_coeffs()
degree = len(coeffs) - 1
for i, coeff in enumerate(coeffs):
    power = degree - i
    if coeff == 0:
        continue
    # Format the sign
    sign = ""
    if i > 0:
        sign = "+ " if coeff > 0 else "- "
    
    # Get the absolute value of the coefficient
    abs_coeff = abs(coeff)
    
    # Format the coefficient and variable part
    if abs_coeff == 1 and power > 0:
        coeff_str = ""
    else:
        coeff_str = str(abs_coeff)

    if power > 0:
        var_str = f"k^{power}" if power > 1 else "k"
        if coeff_str != "":
            term = f"{sign}{coeff_str} * {var_str}"
        else:
            term = f"{sign}{var_str}"
    else: # constant term
        term = f"{sign}{coeff_str}"
        
    # Correct the sign for the very first term
    if i == 0 and coeff < 0:
        term = "-" + term

    terms.append(term)
    
# Join the terms, removing leading '+' signs from the first term if present
final_equation = " ".join(terms).strip()
if final_equation.startswith("+"):
    final_equation = final_equation[2:]

print(f"P(k) = {final_equation}")