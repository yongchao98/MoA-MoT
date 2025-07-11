# Set the desired integer value for n. You can change this value.
n = 5

print(f"This script calculates the q-difference quotient nabla_q(T^n) for n = {n}.")
print("The general formula is: nabla_q(T^n) = [n]_q * T^(n-1)")
print("where [n]_q is the q-integer, defined as (q^n - 1) / (q - 1).\n")

# Step 1: Calculate the exponent of T
exponent = n - 1

# Step 2: Determine the string representation of the coefficient [n]_q
if n == 0:
    # For n=0, [0]_q = 0
    coefficient_str = "0"
elif n > 0:
    # For n > 0, [n]_q = 1 + q + q^2 + ... + q^(n-1)
    if n == 1:
        # [1]_q = 1
        coefficient_str = "1"
    else:
        # Build the sum string: 1 + q + q^2 + ...
        terms = ["1", "q"] + [f"q^{i}" for i in range(2, n)]
        coefficient_str = f"({ ' + '.join(terms) })"
elif n < 0:
    # For n < 0, let n = -m where m > 0.
    # Then [-m]_q = -(q^-1 + q^-2 + ... + q^-m)
    m = -n
    terms = [f"q^{i}" for i in range(-1, n - 1, -1)]
    term_sum_str = " + ".join(terms)
    if m == 1:
        # Avoids parentheses for a single term, e.g., -q^-1
        coefficient_str = f"-{term_sum_str}"
    else:
        coefficient_str = f"-({term_sum_str})"

# Step 3: Assemble the final expression string
if n == 0:
    # If the coefficient is 0, the whole expression is 0
    final_expression = "0"
else:
    # Format the T part of the expression
    if exponent == 0:
        t_part = ""  # T^0 = 1, so it's omitted
    elif exponent == 1:
        t_part = "T"
    else:
        t_part = f"T^{exponent}"
    
    # Combine coefficient and T part for the final result
    if coefficient_str == "1":
        final_expression = t_part if t_part else "1"
    elif t_part == "":
        final_expression = coefficient_str
    else:
        final_expression = f"{coefficient_str} * {t_part}"

# Print the components and the final resulting equation
print("-----------------------------------------")
print(f"For n = {n}:")
print(f"  - The coefficient [n]_q is: {coefficient_str}")
print(f"  - The exponent of T is n-1: {exponent}")
print("-----------------------------------------")
print(f"\nFinal Result:")
print(f"nabla_q(T^{n}) = {final_expression}")
