import math

def predict_inversion_barrier():
    """
    Predicts the inversion barrier for a PAH molecule based on a series.

    The problem provides data for a series of bowl-shaped polycyclic aromatic hydrocarbons (PAHs):
    1. A molecule with 1 bowl unit has an inversion barrier of 10 kcal/mol.
    2. A molecule with 2 bowl units has an inversion barrier of 49 kcal/mol.
    3. The task is to predict the barrier for a molecule with 3 bowl units.

    Let n be the number of bowl units. We look for a relationship B(n) between n and the barrier B.
    Let's assume a quadratic relationship of the form: B(n) = a * n^2 + c
    
    Using the given data points:
    For n = 1: a * (1)^2 + c = 10  =>   a + c = 10  (Eq. 1)
    For n = 2: a * (2)^2 + c = 49  =>  4a + c = 49  (Eq. 2)

    Solving this system of linear equations:
    Subtract Eq. 1 from Eq. 2: (4a + c) - (a + c) = 49 - 10 => 3a = 39 => a = 13.
    Substitute a = 13 into Eq. 1: 13 + c = 10 => c = -3.
    
    So the derived formula is B(n) = 13 * n^2 - 3.
    We will use this formula to predict the barrier for n = 3.
    """
    
    # Coefficients of the quadratic equation B(n) = a*n^2 + c
    a = 13
    c = -3
    
    # The molecule for which we want to predict the barrier has n=3 units
    n_predict = 3
    
    # Calculate the barrier
    predicted_barrier = a * (n_predict**2) + c
    
    # Print the explanation and the result
    print("Based on the provided data, we identified a quadratic relationship.")
    print(f"The formula for the inversion barrier B as a function of the number of units (n) is:")
    print(f"B(n) = {a} * n^2 + ({c})")
    print("\nWe now predict the barrier for the third molecule, which has n = 3 units.")
    print("Calculation:")
    print(f"B({n_predict}) = {a} * {n_predict}^2 - 3")
    print(f"B({n_predict}) = {a} * {n_predict**2} - 3")
    print(f"B({n_predict}) = {a * (n_predict**2)} - 3")
    print(f"B({n_predict}) = {predicted_barrier}")

    print(f"\nThe predicted inversion barrier for triacenaphtho[...]triphenylene is {predicted_barrier} kcal/mol.")

predict_inversion_barrier()