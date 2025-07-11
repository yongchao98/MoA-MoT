import numpy as np

def solve_inversion_barrier():
    """
    Calculates the inversion barrier for a PAH with n=3 pentagons
    based on data from homologues with n=1 and n=2 pentagons.
    """
    # Data points: (n = number of pentagons, B = inversion barrier in kcal/mol)
    n1, B1 = 1, 10
    n2, B2 = 2, 49
    n_predict = 3

    # We assume a quadratic model: B(n) = a*n^2 + b*n
    # This creates a system of two linear equations to find a and b:
    # 1. a*(n1**2) + b*n1 = B1
    # 2. a*(n2**2) + b*n2 = B2
    
    # We can write this in matrix form A*x = y, where x = [a, b]
    A = np.array([[n1**2, n1], [n2**2, n2]])
    y = np.array([B1, B2])
    
    # Solve the system for a and b
    try:
        a, b = np.linalg.solve(A, y)
    except np.linalg.LinAlgError:
        print("Error: Could not solve the system of equations.")
        return

    # Now, use the derived formula to predict the barrier for n=3
    
    print("The relationship between the number of pentagons (n) and the inversion barrier (B) is modeled by the equation:")
    print(f"B(n) = {a:.1f} * n^2 + {b:.1f} * n\n")
    
    print(f"To find the barrier for the molecule with n = {n_predict}, we substitute n into the equation:")
    
    # Print the equation with the numbers substituted, as requested
    term1_val = a * (n_predict**2)
    term2_val = b * n_predict
    final_result = term1_val + term2_val

    print(f"B({n_predict}) = {a:.1f} * {n_predict}^2 + ({b:.1f}) * {n_predict}")
    print(f"B({n_predict}) = {a:.1f} * {n_predict**2} + ({term2_val:.1f})")
    print(f"B({n_predict}) = {term1_val:.1f} + ({term2_val:.1f})")
    print(f"B({n_predict}) = {final_result:.1f}\n")

    # Final answer, rounded to the nearest integer
    predicted_barrier = int(round(final_result))
    
    print(f"The predicted inversion barrier for triacenaphtho[...]triphenylene is {predicted_barrier} kcal/mol.")

if __name__ == '__main__':
    solve_inversion_barrier()