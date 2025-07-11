def generate_nabla_q_expression(n: int) -> str:
    """
    Calculates and returns the string expression for the q-derivative of T^n.

    The q-derivative is given by nabla_q(T^n) = [n]_q * T^(n-1), where
    [n]_q is the q-number. This function expands [n]_q into its sum representation.

    Args:
        n: An integer.

    Returns:
        A string representing the mathematical expression.
    """
    if not isinstance(n, int):
        raise TypeError("Input 'n' must be an integer.")

    # Case n = 0: nabla_q(T^0) = 0
    if n == 0:
        return "0"

    # Determine the exponent of T
    power = n - 1

    # Build the T part of the expression
    if power == 0:
        t_part = ""
    elif power == 1:
        t_part = " * T"
    else:
        t_part = f" * T^({power})"

    # Build the coefficient part ([n]_q)
    coeff_part = ""
    if n > 0:
        # For n > 0, [n]_q = 1 + q + ... + q^(n-1)
        if n == 1:
            return "1"  # [1]_q = 1, T^0 = 1
        
        terms = []
        for k in range(n):
            if k == 0:
                terms.append("1")
            elif k == 1:
                terms.append("q")
            else:
                terms.append(f"q^{k}")
        coeff_part = f"({' + '.join(terms)})"
    else:  # n < 0
        # For n < 0, [n]_q = -(q^n + ... + q^-1)
        terms = []
        for k in range(n, 0):
            if k == -1:
                terms.append("q^-1")
            else:
                terms.append(f"q^({k})")
        
        if len(terms) == 1:
            coeff_part = f"-{terms[0]}"
        else:
            coeff_part = f"-({' + '.join(terms)})"

    return f"{coeff_part}{t_part}"

def main():
    """
    Main function to get user input and print the expression.
    """
    try:
        n_str = input("Please enter an integer value for n: ")
        n = int(n_str)
        expression = generate_nabla_q_expression(n)
        print(f"\nFor n = {n}, the expression for nabla_q(T^n) is:")
        print(expression)
    except (ValueError, TypeError):
        print("Invalid input. Please enter a valid integer.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()