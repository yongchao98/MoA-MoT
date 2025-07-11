import numpy as np

def solve_inversion_barrier():
    """
    This function predicts the inversion barrier for a large PAH molecule
    based on the trend observed from two smaller, related molecules.
    """
    # Known data points
    # N1: number of units for the first molecule (corannulene)
    # B1: inversion barrier for the first molecule
    N1 = 1
    B1 = 10  # kcal/mol

    # N2: number of units for the second molecule (diacenaphtho[...])
    # B2: inversion barrier for the second molecule
    N2 = 2
    B2 = 49  # kcal/mol

    # We assume a quadratic relationship: B = a * N^2 + b
    # This gives a system of two linear equations:
    # 1. a * N1^2 + b = B1
    # 2. a * N2^2 + b = B2

    # We can solve this system for a and b.
    # From eq 1: b = B1 - a * N1^2
    # Substitute into eq 2: a * N2^2 + (B1 - a * N1^2) = B2
    # a * (N2^2 - N1^2) = B2 - B1
    # a = (B2 - B1) / (N2^2 - N1^2)

    a = (B2 - B1) / (N2**2 - N1**2)
    b = B1 - a * (N1**2)

    # Now we have the formula. Let's predict the barrier for the third molecule.
    # N3: number of units for the third molecule (triacenaphtho[...])
    N3 = 3

    # Calculate the predicted barrier B3
    B3 = a * N3**2 + b

    # Print the explanation and the result
    print("Based on the provided data, a quadratic relationship between the number of units (N)")
    print("and the inversion barrier (B) was found of the form: B = a * N^2 + b.")
    print("\nUsing the data points (N=1, B=10) and (N=2, B=49), the coefficients are:")
    print(f"a = {int(a)}")
    print(f"b = {int(b)}")
    print(f"\nThe derived formula is: B = {int(a)} * N^2 - {abs(int(b))}")

    print("\nFor the target molecule, triacenaphtho[...], N = 3.")
    print("The predicted inversion barrier is calculated as follows:")
    
    # Print the final equation with numbers
    print(f"{int(a)} * {N3}^2 - {abs(int(b))} = {int(a)} * {N3**2} - {abs(int(b))} = {int(B3)}")

    print(f"\nThe predicted inversion barrier is {int(B3)} kcal/mol.")

solve_inversion_barrier()