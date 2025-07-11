import numpy as np

def solve_inversion_barrier():
    """
    Calculates the inversion barrier for a target molecule based on a series of related molecules.
    """
    # --- Step 1: Define the known data from the problem ---
    # Let 'n' be the number of "acenaphtho" units, based on the molecule's name and structure.
    # Let 'B' be the corresponding inversion barrier in kcal/mol.

    # Data point 1: dibenzo[ghi,mno]fluoranthene
    n1 = 1
    B1 = 10

    # Data point 2: diacenaphtho[...]chrysene
    n2 = 2
    B2 = 49

    # --- Step 2: Set up and solve the system of equations for the quadratic model B(n) = a*n^2 + b*n ---
    # Equation 1: a*(1)^2 + b*(1) = 10  =>   a +  b = 10
    # Equation 2: a*(2)^2 + b*(2) = 49  =>  4a + 2b = 49
    
    # We represent this system in matrix form A*x = y and solve for x = [a, b]
    A = np.array([[n1**2, n1], [n2**2, n2]])
    y = np.array([B1, B2])
    
    # Solve for the coefficients [a, b]
    try:
        coefficients = np.linalg.solve(A, y)
        a = coefficients[0]
        b = coefficients[1]
    except np.linalg.LinAlgError:
        print("The system of equations could not be solved. The data may be inconsistent.")
        return

    # --- Step 3: Predict the barrier for the target molecule ---
    # The target molecule is triacenaphtho[...]triphenylene, so we set n = 3.
    n3 = 3
    
    # Calculate the predicted barrier B(3) using the derived formula
    B3 = a * n3**2 + b * n3

    # --- Step 4: Print the comprehensive results ---
    print("Based on the provided data, a quadratic model can be used to predict the inversion barrier.")
    print("The model assumes the inversion barrier B(n) is a function of 'n', the number of acenaphtho units.")
    print(f"The relationship is modeled by the formula: B(n) = a*n^2 + b*n")
    print("\nUsing the data points (n=1, B=10) and (n=2, B=49), we solve for the coefficients.")
    print(f"The calculated coefficients are: a = {a}, b = {b}")
    
    print(f"\nThus, the predictive formula is: B(n) = {a}*n^2 + ({b})*n")
    
    print(f"\nFor the target molecule, triacenaphtho[...]triphenylene, n = {n3}.")
    print("The predicted inversion barrier is calculated as follows:")
    
    # We display each number in the final equation as requested
    term1_val = a * n3**2
    term2_val = b * n3
    print(f"B({n3}) = {a} * ({n3}^2) + ({b}) * {n3}")
    print(f"B({n3}) = {a} * {n3**2} + ({term2_val})")
    print(f"B({n3}) = {term1_val} - {abs(term2_val)}")
    
    # The final answer is requested with an accuracy of integers
    final_answer = int(round(B3))
    print(f"B({n3}) = {final_answer}")
    
    print(f"\nThe predicted inversion barrier for triacenaphtho[...]triphenylene is {final_answer} kcal/mol.")


solve_inversion_barrier()