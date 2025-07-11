import sys

def solve_inversion_barrier():
    """
    This script predicts the inversion barrier for the third molecule in a series.

    The problem is modeled by finding a quadratic relationship B(n) = a*n^2 + c that
    fits the data for the first two molecules and using it to predict the value for the third.
    
    Given data:
    - Molecule 1 (n=1): B(1) = 10 kcal/mol
    - Molecule 2 (n=2): B(2) = 49 kcal/mol
    
    We need to find B(3).

    System of equations:
    1) a * (1^2) + c = 10  =>   a + c = 10
    2) a * (2^2) + c = 49  =>  4a + c = 49

    Solving the system:
    Subtract (1) from (2): (4a + c) - (a + c) = 49 - 10 => 3a = 39 => a = 13.
    Substitute a=13 into (1): 13 + c = 10 => c = -3.
    
    The derived formula is B(n) = 13*n^2 - 3.
    This script will verify the formula and calculate B(3).
    """

    # Coefficients of the quadratic relationship B(n) = a*n^2 + c
    a = 13
    c = -3

    # Define the molecule indices and their known barriers
    n1, b1 = 1, 10
    n2, b2 = 2, 49
    n3 = 3

    # Verify the formula with the given data
    b1_calculated = a * n1**2 + c
    b2_calculated = a * n2**2 + c
    
    print("The relationship between molecule index (n) and the inversion barrier B(n) can be described by B(n) = 13*n^2 - 3.")
    print("\nVerification of the formula with the provided data:")
    print(f"For molecule 1 (n={n1}): 13 * {n1}^2 - 3 = {b1_calculated} kcal/mol (matches the given {b1} kcal/mol)")
    print(f"For molecule 2 (n={n2}): 13 * {n2}^2 - 3 = {b2_calculated} kcal/mol (matches the given {b2} kcal/mol)")

    # Calculate the barrier for the third molecule
    b3_calculated = a * n3**2 + c
    
    print("\nPrediction for the third molecule (n=3) using the established formula:")
    # Print the full calculation for the final equation as requested
    print(f"B(3) = 13 * {n3}^2 - 3 = 13 * {n3**2} - 3 = {a * n3**2} - 3 = {b3_calculated}")

    print(f"\nThe predicted inversion barrier is {b3_calculated} kcal/mol.")


solve_inversion_barrier()