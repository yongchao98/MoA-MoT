import numpy as np

def solve_and_predict_barrier():
    """
    Solves for the inversion barrier of the third molecule based on the first two.
    """
    # Step 1: Define the known data points
    # n: number of structural units
    # B: inversion barrier in kcal/mol
    data_points = {
        1: 10,  # n=1, B=10
        2: 49   # n=2, B=49
    }
    
    n_values = list(data_points.keys())
    B_values = list(data_points.values())

    # Step 2: Set up the system of linear equations for B(n) = a*n^2 + b*n
    # a*n1^2 + b*n1 = B1
    # a*n2^2 + b*n2 = B2
    A = np.array([
        [n_values[0]**2, n_values[0]],
        [n_values[1]**2, n_values[1]]
    ])
    
    B = np.array(B_values)

    # Step 3: Solve for the coefficients a and b
    try:
        coeffs = np.linalg.solve(A, B)
        a, b = coeffs[0], coeffs[1]
    except np.linalg.LinAlgError:
        print("Could not solve the system of equations.")
        return

    # The number of units for the molecule in question
    n_predict = 3
    
    # Step 4: Calculate intermediate values and the final prediction for n=3
    n_predict_sq = n_predict**2
    term1 = a * n_predict_sq
    term2 = b * n_predict
    predicted_barrier = int(round(term1 + term2))

    # --- Output the results ---
    print("Based on the provided data, we can model the inversion barrier B as a function of the number of structural units, n.")
    print("Assuming a quadratic relationship B(n) = a*n^2 + b*n:")
    print(f"\nThe derived equation is: B(n) = {a:.1f}*n^2 + ({b:.1f})*n")
    
    print("\nTo find the inversion barrier for the third compound (n=3), we calculate:")
    print(f"B({n_predict}) = {a:.1f} * {n_predict}^2 + ({b:.1f}) * {n_predict}")
    print(f"     = {a:.1f} * {n_predict_sq} - {abs(b):.1f} * {n_predict}")
    print(f"     = {term1:.1f} - {abs(term2):.1f}")
    print(f"     = {predicted_barrier}")
    print(f"\nThe predicted inversion barrier for triacenaphtho[...]triphenylene is {predicted_barrier} kcal/mol.")

if __name__ == "__main__":
    solve_and_predict_barrier()