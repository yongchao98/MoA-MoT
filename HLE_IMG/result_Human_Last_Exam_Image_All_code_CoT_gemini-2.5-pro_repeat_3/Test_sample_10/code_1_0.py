import numpy as np

def solve_and_predict_barrier():
    """
    Solves for the parameters of a quadratic model based on two data points
    and predicts the value for a third point.
    """
    # Data points from the problem description
    # n: index of the molecule in the series
    # b: inversion barrier in kcal/mol
    n1, b1 = 1, 10
    n2, b2 = 2, 49
    n3 = 3

    # We assume a quadratic relationship: B(n) = a*n^2 + c
    # This gives us a system of two linear equations:
    # a*(n1^2) + c = b1  => 1a + c = 10
    # a*(n2^2) + c = b2  => 4a + c = 49

    # We can solve this system using linear algebra.
    # Define the matrix A and vector B for the system Ax = B, where x = [a, c]
    A = np.array([[n1**2, 1], [n2**2, 1]])
    B = np.array([b1, b2])

    # Solve for x = [a, c]
    try:
        a, c = np.linalg.solve(A, B)
    except np.linalg.LinAlgError:
        print("Could not solve the system of equations.")
        return

    # Round coefficients to nearest integers if they are very close
    a = round(a)
    c = round(c)

    # Now, calculate the predicted barrier for molecule 3 (n=3)
    b3 = a * (n3**2) + c

    print("Step 1: Define the data points from the problem.")
    print(f"Molecule 1 (n=1): Inversion Barrier = {b1} kcal/mol")
    print(f"Molecule 2 (n=2): Inversion Barrier = {b2} kcal/mol")
    print("\nStep 2: Establish a mathematical model.")
    print("Assuming a quadratic relationship B(n) = a*n^2 + c, we solve for 'a' and 'c'.")
    print(f"The calculated coefficients are: a = {a}, c = {c}")
    print(f"So, the model is: B(n) = {a}*n^2 + {c}")

    print("\nStep 3: Predict the inversion barrier for Molecule 3 (n=3).")
    print("Plugging n=3 into the equation:")
    print(f"Barrier = {a} * ({n3})^2 + ({c})")
    print(f"Barrier = {a} * ({n3**2}) + ({c})")
    print(f"Barrier = {a * n3**2} + ({c})")
    print(f"Barrier = {b3}")

    print(f"\nThe predicted inversion barrier for triacenaphtho[...]triphenylene is {b3} kcal/mol.")


solve_and_predict_barrier()