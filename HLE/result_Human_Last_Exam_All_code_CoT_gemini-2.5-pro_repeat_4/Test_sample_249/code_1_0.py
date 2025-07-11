import sys

def solve():
    """
    Calculates the minimum possible diameter of a tree with n+2 vertices and m leaves.
    """
    try:
        # Read n and m from command line arguments or use default examples.
        if len(sys.argv) > 2:
            n = int(sys.argv[1])
            m = int(sys.argv[2])
        else:
            # You can change these default values for n and m.
            print("Usage: python your_script.py <n> <m>")
            print("Using example values: n=8, m=3")
            n = 8
            m = 3
    except (ValueError, IndexError):
        print("Invalid input. Please provide two positive integers for n and m.")
        return

    if n <= 0 or m <= 0:
        print("n and m must be positive integers.")
        return

    # A tree with n+2 vertices can have at most n+1 leaves.
    if m > n + 1:
        print(f"For n={n}, m={m}:")
        print("A tree with n+2 vertices cannot have more than n+1 leaves. This case is impossible.")
        return
    
    print(f"For n={n}, m={m}:")

    # Case 1: m = n + 1 (implies 1 internal vertex) -> Star graph
    if m == n + 1:
        diameter = 2
        print("The minimum possible diameter is 2.")

    # Case 2: m = n (implies 2 internal vertices) -> Dumbbell shape
    elif m == n:
        diameter = 3
        print("The minimum possible diameter is 3.")

    # Case 3: 2m < n + 1 (leaf-starved) -> Path-like internal subgraph
    elif 2 * m < n + 1:
        diameter = n - m + 3
        print(f"The minimum possible diameter is n - m + 3 = {n} - {m} + 3 = {diameter}.")

    # Case 4: m <= n-1 and 2m >= n+1 (leaf-rich) -> Star-like internal subgraph
    else:
        diameter = 4
        print("The minimum possible diameter is 4.")

if __name__ == '__main__':
    solve()
