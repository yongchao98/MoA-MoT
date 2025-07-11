import sys

def solve_q_derivative_of_T_n(n: int):
    """
    This function computes and prints the symbolic formula for the q-derivative
    of T^n. It also outputs the integer numbers that define the structure of
    the resulting equation, as requested.

    Args:
        n: A non-negative integer exponent.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: The exponent n must be a non-negative integer.", file=sys.stderr)
        return

    # --- Construct the symbolic equation ---

    # Left-hand side of the equation
    lhs = f"Nabla_q(T^{n})"

    # Right-hand side of the equation
    if n == 0:
        # For n=0, Nabla_q(T^0) = Nabla_q(1) = (1-1)/((q-1)T) = 0
        rhs = "0"
    elif n == 1:
        # For n=1, Nabla_q(T^1) = [1]_q * T^0 = 1 * 1 = 1
        rhs = "1"
    else:
        # For n > 1, construct the q-polynomial for [n]_q
        q_terms = []
        for i in range(n):
            if i == 0:
                q_terms.append("1")
            elif i == 1:
                q_terms.append("q")
            else:
                q_terms.append(f"q^{i}")
        q_poly_str = "(" + " + ".join(q_terms) + ")"

        # Construct the T part, T^(n-1)
        t_exponent = n - 1
        if t_exponent == 1:
            t_part_str = "T"
        else:
            t_part_str = f"T^{t_exponent}"
        
        rhs = f"{q_poly_str} * {t_part_str}"
        
    print("The formula for the q-derivative of T^n is:")
    print(f"{lhs} = {rhs}\n")

    # --- Output the numbers in the final equation ---
    
    print("The numbers that define the final equation are:")
    
    # Exponent of T on the left-hand side
    print(f"LHS T exponent: {n}")

    if n == 0:
        print("RHS value: 0")
    else:
        # [n]_q polynomial coefficients (all are 1)
        q_poly_coeffs = [1] * n
        print(f"q-poly coefficients: {' '.join(map(str, q_poly_coeffs))}")
        
        # [n]_q polynomial exponents
        q_poly_exponents = list(range(n))
        print(f"q-poly exponents: {' '.join(map(str, q_poly_exponents))}")

        # Exponent of T on the right-hand side
        rhs_t_exponent = n - 1
        print(f"RHS T exponent: {rhs_t_exponent}")


# --- Example execution with n = 5 ---
if __name__ == '__main__':
    # You can change this value to see the result for different n
    n_example = 5
    print(f"--- Calculating for n = {n_example} ---")
    solve_q_derivative_of_T_n(n_example)