import sys

def solve_nabla_q(n):
    """
    This function computes and prints the symbolic expression for the q-difference
    quotient of T^n, denoted as nabla_q(T^n).

    The calculation is based on the formula:
    nabla_q(T^n) = [n]_q * T^(n-1)
    where [n]_q is the q-number of n.
    """
    if not isinstance(n, int):
        print("Error: Input 'n' must be an integer.", file=sys.stderr)
        return

    # Format the left-hand side of the equation
    equation_lhs = f"nabla_q(T^{n})"
    if n < 0:
        # Add parentheses for negative n for clarity
        equation_lhs = f"nabla_q(T^({n}))"

    # Case 1: n = 0
    # nabla_q(T^0) = nabla_q(1) = 0
    if n == 0:
        print(f"{equation_lhs} = 0")
        return

    # Case 2: n = 1
    # nabla_q(T^1) = [1]_q * T^0 = 1 * 1 = 1
    if n == 1:
        print(f"{equation_lhs} = 1")
        return
        
    # Determine the T part of the expression
    exponent_T = n - 1
    if exponent_T == 1:
        t_part = " * T"
    else:
        # Use parentheses for exponents other than 1
        t_part = f" * T^({exponent_T})"

    # Case 3: n > 1
    if n > 1:
        # [n]_q = 1 + q + q^2 + ... + q^(n-1)
        q_terms = []
        for i in range(n):
            if i == 0:
                q_terms.append("1")
            elif i == 1:
                q_terms.append("q")
            else:
                q_terms.append(f"q^{i}")
        q_poly_str = " + ".join(q_terms)
        
        # Final expression for n > 1
        print(f"{equation_lhs} = ({q_poly_str}){t_part}")

    # Case 4: n < 0
    elif n < 0:
        # [n]_q = -q^n * [|n|]_q
        m = -n  # m is a positive integer
        q_terms = []
        for i in range(m):
            if i == 0:
                q_terms.append("1")
            elif i == 1:
                q_terms.append("q")
            else:
                q_terms.append(f"q^{i}")
        q_poly_str = " + ".join(q_terms)
        
        # Add parentheses if the polynomial has more than one term
        if m > 1:
            q_poly_str = f"({q_poly_str})"
        
        # Final expression for n < 0
        print(f"{equation_lhs} = -q^({n}) * {q_poly_str}{t_part}")


def main():
    """
    Main function to parse command-line arguments and run the solver.
    """
    # Set a default value for n
    n_default = 4
    n = n_default

    if len(sys.argv) > 1:
        try:
            # Use the integer provided by the user
            n = int(sys.argv[1])
        except (ValueError, IndexError):
            print(f"Usage: python {sys.argv[0]} [integer n]", file=sys.stderr)
            print(f"Invalid input. Using default value n={n_default}.", file=sys.stderr)
            n = n_default

    print(f"Calculating nabla_q(T^n) for n = {n}:\n")
    solve_nabla_q(n)

if __name__ == "__main__":
    main()