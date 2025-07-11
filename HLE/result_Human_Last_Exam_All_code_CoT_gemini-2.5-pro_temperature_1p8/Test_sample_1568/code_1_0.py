# Define constants from the problem statement.
power = 3
n_start = 3
excluded_indices = [1, 2]

# Define Unicode symbols for clearer mathematical notation in the output.
prod_sym = "\u220F"      # Product symbol (Π)
inf_sym = "\u221E"       # Infinity symbol (∞)
gamma_sym = "\u0393"     # Gamma symbol (Γ)
omega_sym = "\u03C9"     # Omega symbol (ω)
pi_sym = "\u03C0"        # Pi symbol (π)

# Construct the left-hand side of the equation string.
lhs = f"{prod_sym}_{{n={n_start}}}^{inf_sym} (1 - z^{power}/n^{power})"

# Construct the right-hand side of the equation string.

# First, create the part for the terms that are divided out (n=1 and n=2).
excluded_terms_list = []
for n in excluded_indices:
    denominator = n**power
    if denominator == 1:
        excluded_terms_list.append(f"(1 - z^{power})")
    else:
        excluded_terms_list.append(f"(1 - z^{power}/{denominator})")
excluded_part_str = " * ".join(excluded_terms_list)

# Then, create the Gamma function product part.
gamma_part_str = f"{gamma_sym}(1-z){gamma_sym}(1 - z*{omega_sym}){gamma_sym}(1 - z*{omega_sym}^2)"

# Assemble the full right-hand side.
rhs = f"1 / ({excluded_part_str} * {gamma_part_str})"

# Combine LHS and RHS to form the final equation.
final_equation = f"{lhs} = {rhs}"

# Create the definition for the omega symbol.
omega_definition = f"where {omega_sym} = exp(2*{pi_sym}*i/{power}) is a primitive cube root of unity."

# Print the final result.
print(final_equation)
print(omega_definition)