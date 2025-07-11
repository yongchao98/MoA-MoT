def print_q_nabla_formula(n: int):
    """
    This function calculates and prints the formula for nabla_q(T^n)
    for a given non-negative integer n.
    
    Args:
        n (int): A non-negative integer representing the power of T.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: Please provide a non-negative integer for n.")
        return

    # The left-hand side of the equation
    lhs = f"nabla_q(T^{n})"

    # Handle the special cases for n=0 and n=1
    if n == 0:
        # nabla_q(T^0) = nabla_q(1) = (1-1)/(q-1)T = 0
        rhs = "0"
    elif n == 1:
        # nabla_q(T^1) = [1]_q * T^0 = 1 * 1 = 1
        rhs = "1"
    else:
        # Build the coefficient part: [n]_q = (1 + q + q^2 + ... + q^(n-1))
        q_terms = []
        for i in range(n):
            if i == 0:
                q_terms.append("1")
            elif i == 1:
                q_terms.append("q")
            else:
                q_terms.append(f"q^{i}")
        
        coeff_str = f"({ ' + '.join(q_terms) })"

        # Build the T part: T^(n-1)
        exponent = n - 1
        if exponent == 1:
            t_part = "T"
        else:
            t_part = f"T^{exponent}"
        
        rhs = f"{coeff_str} * {t_part}"
    
    # Print the final, formatted equation
    print(f"{lhs} = {rhs}")

# --- Main execution ---
# Set the integer n for which to generate the formula.
# You can change this value to see the formula for different powers.
n_value = 5
print_q_nabla_formula(n_value)