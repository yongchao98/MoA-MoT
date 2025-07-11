def solve():
    """
    Calculates the minimum possible diameter of a tree with n+2 vertices and m leaves.
    """
    print("This program calculates the minimum possible diameter of a tree with n+2 vertices and m leaves.")
    try:
        n_str = input("Enter the value for n (a positive integer): ")
        m_str = input("Enter the value for m (a positive integer): ")
        n = int(n_str)
        m = int(m_str)

        if n <= 0 or m <= 0:
            raise ValueError("n and m must be positive.")

        # A tree with V=n+2 vertices has at least 2 leaves (if V>=3) and at most V-1=n+1 leaves.
        if n + 2 >= 3 and m < 2:
            print(f"Error: A tree with {n+2} vertices must have at least 2 leaves, but m={m}.")
            return
        if m > n + 1:
            print(f"Error: A tree with {n+2} vertices has at most {n+1} leaves, but m={m}.")
            return
        # If m=n, we need n>=2 (implies m>=2) for a valid tree structure described in the derivation.
        # If n=1, V=3, the only possible tree is P_3 with m=2=n+1.
        if m == n and n < 2:
             print(f"Error: If n=1, m must be 2, not 1.")
             return

        # Calculate and print the diameter based on the case analysis
        if m == n + 1:
            diameter = 2
            print("Case: m = n + 1")
            print(f"The minimum possible diameter is 2.")
        elif m == n:
            diameter = 3
            print("Case: m = n")
            print(f"The minimum possible diameter is 3.")
        elif m < n:
            if 2 * m >= n + 1:
                diameter = 4
                print("Case: m < n and 2*m >= n+1")
                print(f"The minimum possible diameter is 4.")
            else:  # 2 * m < n + 1
                diameter = n - 2 * m + 5
                print("Case: m < n and 2*m < n+1")
                print(f"The minimum possible diameter is calculated by the formula: n - 2*m + 5")
                print(f"Final Equation: {n} - 2 * {m} + 5 = {diameter}")

    except ValueError as e:
        print(f"Invalid input. Please enter positive integers for n and m. Details: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

solve()