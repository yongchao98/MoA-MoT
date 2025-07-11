import sys

def solve():
    """
    Calculates the minimum possible diameter of an undirected tree with n+2 vertices and m leaves.
    """
    try:
        # Prompt the user for n and m
        n_str = input("Enter the positive integer n: ")
        m_str = input("Enter the positive integer m: ")
        
        n = int(n_str)
        m = int(m_str)
        
        if n <= 0 or m <= 0:
            print("Error: n and m must be positive integers.")
            return

        # A tree needs at least 2 leaves, unless it's a single vertex.
        # For a diameter to exist, we need at least 2 vertices.
        if m < 2:
            print("Error: A tree with a diameter must have at least 2 leaves (m >= 2).")
            return

        print("-" * 20)
        
        # Case 1: n <= m (few internal nodes)
        if n <= m:
            result = n + 3 - m
            print(f"Since n <= m ({n} <= {m}), the minimum diameter is given by the formula n + 3 - m.")
            print(f"The calculation is: {n} + 3 - {m} = {result}")
            final_answer = result
    
        # Case 2: n > m (many internal nodes)
        else:
            # Subcase 2a: "leaf-rich"
            if 2 * m >= n + 1:
                result = 4
                print(f"Since n > m ({n} > {m}) and 2*m >= n+1 ({2*m} >= {n+1}), there are relatively many leaves.")
                print(f"The minimum possible diameter is 4.")
                final_answer = result
            # Subcase 2b: "leaf-poor"
            else:
                result = n - 2 * m + 5
                print(f"Since n > m ({n} > {m}) and 2*m < n+1 ({2*m} < {n+1}), there are relatively few leaves.")
                print(f"The minimum diameter is given by the formula n - 2*m + 5.")
                print(f"The calculation is: {n} - 2 * {m} + 5 = {result}")
                final_answer = result

        print("-" * 20)
        # We wrap the final answer as requested
        # print(f"<<<{final_answer}>>>")
        
    except ValueError:
        print("Invalid input. Please enter valid integers for n and m.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)

solve()