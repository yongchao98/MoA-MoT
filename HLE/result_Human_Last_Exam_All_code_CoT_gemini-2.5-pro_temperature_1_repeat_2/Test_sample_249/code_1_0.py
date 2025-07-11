import math

def solve():
    """
    Calculates the minimum possible diameter of a tree with n+2 vertices and m leaves.
    """
    try:
        # In a real scenario, n and m would be inputs.
        # For this example, let's use some placeholder values.
        # The user can replace these with the desired n and m.
        # To make the code executable, we'll hardcode some values.
        # Let's use n=6, m=3 as an example from the thought process.
        n_str = "6"
        m_str = "3"
        # In a real use case, you might get input like this:
        # n_str = input("Enter n: ")
        # m_str = input("Enter m: ")
        
        n = int(n_str)
        m = int(m_str)

        if n <= 0 or m <= 0:
            print("n and m must be positive integers.")
            return

        # A tree with V > 2 vertices must have at least 2 leaves.
        # A tree with V = n+2 vertices can have at most n+1 leaves (a star graph).
        # The number of internal nodes must be non-negative (n+2 >= m).
        if m < 2 and n + 2 > 2:
            print(f"For a tree with {n+2} vertices (which is > 2), m must be >= 2.")
            return
        if m > n + 1:
            print(f"For a tree with {n+2} vertices, m cannot exceed {n+1}.")
            return


        # Case 1: 2*m is greater than or equal to n+1
        if 2 * m >= n + 1:
            if m == n + 1:
                diameter = 2
                print(f"Given n={n} and m={m}:")
                print("The condition 2*m >= n+1 is met.")
                print("Since m = n+1, the tree is a star graph.")
                print(f"The minimum diameter is 2.")
            elif m == n:
                diameter = 3
                print(f"Given n={n} and m={m}:")
                print("The condition 2*m >= n+1 is met.")
                print("Since m = n, the minimal structure is a bi-star.")
                print(f"The minimum diameter is 3.")
            else: # 2*m >= n+1 and m <= n-1
                diameter = 4
                print(f"Given n={n} and m={m}:")
                print("The condition 2*m >= n+1 is met and m < n.")
                print("The minimal structure has a star-like core of internal nodes.")
                print(f"The minimum diameter is 4.")
        # Case 2: 2*m is less than n+1
        else:
            diameter = n + 5 - 2 * m
            print(f"Given n={n} and m={m}:")
            print("The condition 2*m < n+1 is met.")
            print("The tree is 'leaf-poor', forcing a more path-like core.")
            print(f"The minimum diameter is calculated by the formula: n + 5 - 2*m")
            print(f"Result: {n} + 5 - 2*{m} = {n + 5 - 2 * m}")

    except (ValueError, TypeError):
        print("Please provide valid positive integers for n and m.")

solve()