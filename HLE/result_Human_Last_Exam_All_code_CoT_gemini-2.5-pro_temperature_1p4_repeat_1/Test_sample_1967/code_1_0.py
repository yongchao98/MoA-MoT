def solve_destabilizers():
    """
    Calculates the number of different destabilizer sets for the n-qubit
    stabilizer generator set {Z_1, ..., Z_n}.
    """
    try:
        n_str = input("Enter the number of qubits (n): ")
        n = int(n_str)
        if n < 1:
            print("Error: The number of qubits (n) must be a positive integer.")
            return

        # The formula for the number of sets is 2^((n^2 + 5n) / 2).
        # We calculate the components of the formula step-by-step as requested.
        
        n_squared = n * n
        five_n = 5 * n
        numerator = n_squared + five_n
        
        # The numerator n*(n+5) is always even, so the exponent is an integer.
        exponent = numerator // 2
        
        # Python's integers handle arbitrary size, so this won't overflow.
        result = 2**exponent

        print("\n--- Calculation Steps ---")
        print(f"For n = {n}, the formula for the number of sets is 2^((n^2 + 5n) / 2).")
        print(f"1. Calculate the terms in the exponent's numerator:")
        print(f"   n^2 = {n}^2 = {n_squared}")
        print(f"   5n = 5 * {n} = {five_n}")
        print(f"2. Sum them up:")
        print(f"   n^2 + 5n = {n_squared} + {five_n} = {numerator}")
        print(f"3. Divide by 2 to get the final exponent:")
        print(f"   ({numerator}) / 2 = {exponent}")
        print(f"4. Calculate the final result:")
        print(f"   2^{exponent}")

        print("\n--- Final Answer ---")
        print(f"The total number of different destabilizer sets is: {result}")

    except ValueError:
        print("Error: Invalid input. Please enter a valid integer for n.")
    except Exception as e:
        print(f"An error occurred: {e}")

solve_destabilizers()