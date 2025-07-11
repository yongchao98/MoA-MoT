import sys

def solve():
    """
    Calculates the maximum value of n in terms of k for a k-uniform intersecting family
    with full differences of size k-1.
    The formula is n = k^2 - k + 1.
    """
    if len(sys.argv) != 2:
        print("Usage: python your_script.py <k>")
        return

    try:
        k = int(sys.argv[1])
        if k < 2:
            print("Error: k must be an integer greater than or equal to 2.")
            return
    except ValueError:
        print("Error: k must be an integer.")
        return

    # The maximum value of n is given by the formula k^2 - k + 1.
    n = k**2 - k + 1
    
    print(f"For k = {k}, the maximum value of n is given by the expression:")
    print(f"n = k^2 - k + 1")
    print(f"n = {k}^2 - {k} + 1")
    print(f"n = {k*k} - {k} + 1")
    print(f"n = {k*k - k} + 1")
    print(f"n = {n}")


if __name__ == "__main__":
    solve()
