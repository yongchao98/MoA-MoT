import numpy as np

def solve_inversion_barrier():
    """
    Solves for the inversion barrier of the third molecule based on the first two.
    """
    # Step 1 & 2: Define the known data points based on the problem description.
    # n: index of the molecule in the series
    # barriers: corresponding inversion barrier in kcal/mol
    n_known = np.array([1, 2])
    barriers_known = np.array([10, 49])

    # Step 3: Assume a quadratic model B(n) = A*n^2 + B*n.
    # This creates a system of linear equations:
    # For n=1: A*(1)^2 + B*(1) = 10  =>  1*A + 1*B = 10
    # For n=2: A*(2)^2 + B*(2) = 49  =>  4*A + 2*B = 49
    # We can write this in matrix form M * x = v, where x = [A, B].
    
    M = np.array([[n**2, n] for n in n_known])
    v = barriers_known

    # Step 4: Solve the system of equations for the coefficients A and B.
    try:
        coeffs = np.linalg.solve(M, v)
        A = coeffs[0]
        B = coeffs[1]
    except np.linalg.LinAlgError:
        print("Could not solve the system of equations.")
        return

    # Step 5: Calculate the predicted barrier for the target molecule (n=3).
    n_target = 3
    predicted_barrier = A * n_target**2 + B * n_target

    # Print the explanation and the final result.
    print("A quadratic relationship B(n) = A*n^2 + B*n is assumed, where n is the molecule index.")
    print(f"Using the data B(1) = {barriers_known[0]} and B(2) = {barriers_known[1]}, we solve for A and B.")
    print(f"The calculated coefficients are A = {A} and B = {B}.")
    print("\nThe derived equation for the inversion barrier is:")
    print(f"Barrier(n) = {A} * n^2 + {B} * n")
    
    print("\nUsing this equation to predict the barrier for the target molecule (n=3):")
    
    # Show the calculation step-by-step as requested.
    print(f"Barrier({n_target}) = {A} * {n_target}^2 + ({B}) * {n_target}")
    print(f"Barrier({n_target}) = {A} * {n_target**2} + {B * n_target}")
    print(f"Barrier({n_target}) = {A * n_target**2} - {abs(B * n_target)}")
    print(f"Barrier({n_target}) = {int(predicted_barrier)}")

    print(f"\nThe predicted inversion barrier for triacenaphtho[...]triphenylene is {int(predicted_barrier)} kcal/mol.")

solve_inversion_barrier()