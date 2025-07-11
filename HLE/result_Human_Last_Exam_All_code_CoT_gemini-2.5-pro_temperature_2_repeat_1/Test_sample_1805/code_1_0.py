# The integer power n for T^n. You can change this value.
n = 5

if not isinstance(n, int) or n < 0:
    print("Error: n must be a non-negative integer.")
else:
    # Start building the string for the equation
    equation_lhs = f"nabla_q(T^{n})"

    # Handle the cases for n
    if n == 0:
        #nabla_q(T^0) = nabla_q(1) = 0
        equation_rhs = "0"
    elif n == 1:
        #nabla_q(T^1) = (qT - T) / ((q-1)T) = 1
        equation_rhs = "1"
    else: # n > 1
        # Build the coefficient part [n]_q = 1 + q + q^2 + ... + q^(n-1)
        q_terms = ["1", "q"]
        if n > 2:
            q_terms.extend([f"q^{k}" for k in range(2, n)])
        
        # The sum is the coefficient of the T term
        coeff_str = " + ".join(q_terms)
        
        # Build the T part, which is T^(n-1)
        if (n - 1) == 1:
            t_part_str = "T"
        else:
            t_part_str = f"T^{n-1}"
        
        # Combine coefficient and T part
        equation_rhs = f"({coeff_str}){t_part_str}"

    # Print the final equation
    print(f"{equation_lhs} = {equation_rhs}")
