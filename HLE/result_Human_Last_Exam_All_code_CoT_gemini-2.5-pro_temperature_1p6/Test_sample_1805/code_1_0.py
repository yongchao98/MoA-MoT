import sys

def solve_q_derivative(n):
    """
    Calculates and prints the symbolic expression for nabla_q(T^n),
    where nabla_q is the q-difference quotient.
    The final expression is presented as [n]_q * T^(n-1), with the q-analogue [n]_q
    expanded into its polynomial form.
    The function prints the final equation showing all numerical components as requested.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: Input 'n' must be a non-negative integer.", file=sys.stderr)
        return

    # Handle the case n=0, where T^0=1 is a constant and its derivative is 0.
    if n == 0:
        print(f"nabla_q(T^{n}) = 0")
        return

    # Handle the case n=1, where nabla_q(T) = 1.
    if n == 1:
        print(f"nabla_q(T^{n}) = 1")
        return

    # Build the string for the q-analogue [n]_q = 1 + q + q^2 + ... + q^(n-1)
    # Handle n=2 separately for nice formatting (no exponent on the 'q' term).
    if n == 2:
        q_poly = "(1 + q)"
    else:
        # Create a list of terms: "1", "q", "q^2", ...
        terms = ["1", "q"]
        for i in range(2, n):
            terms.append(f"q^{i}")
        q_poly = f"({' + '.join(terms)})"

    # Build the string for the T part, T^(n-1).
    # Handle the exponent being 1 for nice formatting (no "^1").
    if (n - 1) == 1:
        t_poly = "T"
    else:
        t_poly = f"T^{n - 1}"

    # Print the final equation, showing T^n and the resulting expression.
    print(f"nabla_q(T^{n}) = {q_poly} * {t_poly}")


# Set the integer n for the calculation.
# You can change this value to see the result for a different n.
n_integer = 4
solve_q_derivative(n_integer)
