def q_derivative_of_T_n(n: int):
    """
    Calculates and prints the symbolic expression for the q-derivative of T^n.

    The q-derivative is defined as nabla_q(f(T)) = (f(qT) - f(T)) / (qT - T).
    For f(T) = T^n, this results in [n]_q * T^(n-1).

    Args:
        n: An integer.
    """
    if not isinstance(n, int):
        print("Error: n must be an integer.")
        return

    # Case n = 0
    if n == 0:
        # nabla_q(T^0) = nabla_q(1) = 0.
        # This is consistent with [0]_q = 0.
        print(f"nabla_q(T^{n}) = 0")
        return

    # Case n > 0
    if n > 0:
        # The coefficient is the q-number [n]_q = 1 + q + q^2 + ... + q^(n-1)
        q_bracket_terms = []
        for k in range(n):
            if k == 0:
                # Term is 1*q^0 = 1
                q_bracket_terms.append("1")
            elif k == 1:
                # Term is 1*q^1 = q
                q_bracket_terms.append("q")
            else:
                # Term is 1*q^k
                q_bracket_terms.append(f"q^{k}")
        
        q_bracket_str = " + ".join(q_bracket_terms)

        # The power of T is n-1
        power_of_T = n - 1
        if power_of_T == 0:
            # For n=1, T^0 = 1, the result is just [1]_q = 1
            final_expression = q_bracket_str
        elif power_of_T == 1:
            # For n=2, T^1 = T
            final_expression = f"({q_bracket_str}) * T"
        else:
            final_expression = f"({q_bracket_str}) * T^{power_of_T}"
    
    # Case n < 0
    else: # n is negative
        # Let n = -m. The result is [-m]_q * T^(-m-1).
        # [-m]_q = -(q^-m + q^(-m+1) + ... + q^-1)
        
        # The loop runs from n up to -1.
        neg_q_bracket_terms = []
        for k in range(n, 0):
            # Term is q^k
            neg_q_bracket_terms.append(f"q^{k}")
        
        neg_q_bracket_str = " + ".join(neg_q_bracket_terms)

        # The power of T is n-1
        power_of_T = n - 1
        
        final_expression = f"-({neg_q_bracket_str}) * T^{power_of_T}"

    print(f"nabla_q(T^{n}) = {final_expression}")


# --- Demonstration ---
# You can change the value of n to see the result for different integers.
# Example for a positive integer:
print("Demonstration for n = 5:")
q_derivative_of_T_n(5)

# Example for a negative integer:
print("\nDemonstration for n = -4:")
q_derivative_of_T_n(-4)

# Example for n = 1:
print("\nDemonstration for n = 1:")
q_derivative_of_T_n(1)

# Example for n = 0:
print("\nDemonstration for n = 0:")
q_derivative_of_T_n(0)