import sys

def solve():
    """
    Calculates the maximum number of eigenvalues greater than 2 that M(C_n) can have.

    Let n be the number of objects in the category C_n.
    The problem reduces to finding an integer partition of n, n = n_1 + ... + n_k,
    that maximizes the number of certain eigenvalues greater than 2.

    Through analysis of the eigenvalue equations, we find the following:
    - If n = 1, we can't have any eigenvalues greater than 2. The single eigenvalue is 2. Max is 0.
    - If n = 2, partitions are (2) or (1,1).
      - (2): The eigenvalues are n+N = 2+2=4 and 0. One eigenvalue > 2.
      - (1,1): The eigenvalues are n+N = 2+1=3 and 1. One eigenvalue > 2.
      Max is 1.
    - If n >= 3 and odd, let n = 2k+1. The optimal partition is one part of size 3 and k-1 parts of size 2.
      This gives k = (n-1)/2 eigenvalues > 2.
    - If n >= 4 and even, let n = 2k. The optimal partition is two parts of size 3 and k-3 parts of size 2.
      This gives k-1 = n/2 - 1 eigenvalues > 2.

    The code implements this derived formula.
    """
    try:
        n_str = sys.argv[1]
        n = int(n_str)
    except (IndexError, ValueError):
        # Default value for demonstration if no command-line argument is provided.
        n = 8 

    if n <= 0:
        result = 0
    elif n == 1:
        result = 0
    elif n == 2:
        result = 1
    elif n % 2 == 1: # Odd n >= 3
        result = (n - 1) // 2
    else: # Even n >= 4
        result = n // 2 - 1

    print(result)

if __name__ == "__main__":
    solve()
