import numpy as np

def predict_inversion_barrier():
    """
    This script predicts the inversion barrier for a PAH molecule based on a trend
    observed in a series of related molecules.

    The model assumes a quadratic relationship between the molecule index (n)
    and its inversion barrier B(n), of the form B(n) = a*n^2 + b*n.
    The coefficients 'a' and 'b' are determined from the known data points.
    """
    
    # Known data points:
    # n=1: B(1) = 10 kcal/mol
    # n=2: B(2) = 49 kcal/mol
    n_known = np.array([1, 2])
    barriers_known = np.array([10, 49])
    
    # We set up the system of linear equations:
    # a*(1)^2 + b*(1) = 10
    # a*(2)^2 + b*(2) = 49
    # This is equivalent to Ax = B, where x = [a, b]
    A = np.array([[n**2, n] for n in n_known])
    B = barriers_known
    
    # Solve for x = [a, b]
    try:
        a, b = np.linalg.solve(A, B)
    except np.linalg.LinAlgError:
        print("Could not solve the system of equations.")
        return

    # The molecule for which we want to predict the barrier
    n_predict = 3
    
    # Calculate the predicted barrier using the formula B(n) = a*n^2 + b*n
    predicted_barrier = a * n_predict**2 + b * n_predict
    
    print("The relationship between the molecule index 'n' and the inversion barrier 'B(n)' is modeled by the formula:")
    print(f"B(n) = {a:.1f}*n^2 + ({b:.1f})*n\n")
    
    print("To predict the inversion barrier for the target molecule (n=3), we calculate:")
    
    # Print the equation with the final numbers
    # Ensure minus sign is shown correctly for the second term
    if b < 0:
        print(f"{a:.1f} * {n_predict}^2 - {abs(b):.1f} * {n_predict} = {int(predicted_barrier)}")
    else:
        print(f"{a:.1f} * {n_predict}^2 + {b:.1f} * {n_predict} = {int(predicted_barrier)}")
        
    print(f"\nThe predicted inversion barrier is {int(predicted_barrier)} kcal/mol.")

# Execute the function to get the prediction
predict_inversion_barrier()
