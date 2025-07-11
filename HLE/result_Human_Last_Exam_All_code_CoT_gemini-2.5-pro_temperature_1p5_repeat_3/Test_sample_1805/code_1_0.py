def solve_q_nabla(n):
    """
    Generates the symbolic expression for the q-derivative of T^n.
    The formula is: nabla_q(T^n) = [n]_q * T^(n-1).
    This function constructs a string representing this result.
    """
    if not isinstance(n, int):
        return "Error: n must be an integer."

    # Case 1: n = 0
    # The q-derivative of a constant T^0 = 1 is 0.
    if n == 0:
        return "0"

    # Case 2: n > 0 (positive integers)
    if n > 0:
        # Build the q-number part, [n]_q
        if n == 1:
            q_part_str = "1"
        else:
            # For n > 1, [n]_q = 1 + q + q^2 + ... + q^(n-1)
            terms = ["1"]
            if n > 1:
                terms.append("q")
            for k in range(2, n):
                terms.append(f"q^{k}")
            q_part_str = f"({ ' + '.join(terms) })"

        # Build the T part, T^(n-1)
        exponent = n - 1
        if exponent == 0:
            # T^0 = 1, so no T part is needed for n=1
            t_part_str = ""
        elif exponent == 1:
            # T^1 = T
            t_part_str = " * T"
        else:
            t_part_str = f" * T^{exponent}"
        
        return q_part_str + t_part_str
        
    # Case 3: n < 0 (negative integers)
    else: # n is a negative integer
        # The q-number is [n]_q = (q^n - 1) / (q - 1)
        q_part_str = f"((q^({n}) - 1)/(q - 1))"
        
        # The T part is T^(n-1)
        exponent = n - 1
        t_part_str = f" * T^({exponent})"

        return q_part_str + t_part_str


def main():
    """
    Main function to demonstrate the calculation for a specific value of n.
    """
    print("This script calculates the q-derivative of T^n, denoted as nabla_q(T^n).")
    print("The general formula is: nabla_q(T^n) = [n]_q * T^(n-1)\n")
    
    # You can change this integer value to see the result for different powers.
    n = 5 
    
    # Get the symbolic expression from the solver function
    result = solve_q_nabla(n)
    
    # Print the final result in an equation format
    print(f"For n = {n}:")
    
    # Create the left-hand side of the equation string
    if n == 1:
      lhs = "nabla_q(T)"
    else:
      lhs = f"nabla_q(T^{n})"
    
    print(f"{lhs} = {result}")

if __name__ == "__main__":
    main()