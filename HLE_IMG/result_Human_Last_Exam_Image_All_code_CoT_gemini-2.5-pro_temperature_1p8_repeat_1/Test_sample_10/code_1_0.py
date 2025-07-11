import numpy as np

def solve_inversion_barrier():
    """
    This function predicts the inversion barrier of a PAH based on a series of related molecules.
    """
    # Let N be the number of "bowl" units.
    # Data from the problem:
    # N1=1, Barrier1=10
    # N2=2, Barrier2=49
    # N3=3, Barrier3=?

    N1, B1 = 1, 10
    N2, B2 = 2, 49
    N3 = 3

    print("Step 1: Define the relationship between the number of bowl units (N) and the inversion barrier B(N).")
    print("We assume a quadratic model: B(N) = a * N^2 + b * N\n")

    print("Step 2: Use the provided data to set up a system of linear equations.")
    print(f"For N={N1}, B({N1})={B1}  =>  a*({N1}^2) + b*({N1}) = {B1}")
    print(f"For N={N2}, B({N2})={B2}  =>  a*({N2}^2) + b*({N2}) = {B2}\n")

    # System of equations:
    # a + b = 10
    # 4a + 2b = 49
    # We can solve this using numpy's linear algebra solver.
    # Matrix A represents the coefficients of a and b. Vector y represents the results.
    A = np.array([[N1**2, N1], [N2**2, N2]])
    y = np.array([B1, B2])

    # Solve Ax = y for x, where x = [a, b]
    try:
        coefficients = np.linalg.solve(A, y)
        a = coefficients[0]
        b = coefficients[1]
    except np.linalg.LinAlgError:
        print("Could not solve the system of equations.")
        return

    print("Step 3: Solve for the coefficients 'a' and 'b'.")
    print(f"The calculated coefficients are: a = {a}, b = {b}\n")
    
    print("Step 4: Predict the barrier for the third molecule (N=3) using the derived formula.")
    predicted_barrier = a * N3**2 + b * N3
    final_answer = int(round(predicted_barrier))

    # The final output prints the full equation with the numbers plugged in
    print("The final calculation is:")
    print(f"{a} * {N3}^2 + ({b}) * {N3} = {a * N3**2} - {abs(b * N3)} = {final_answer}")
    print(f"\nThe predicted inversion barrier for triacenaphtho[...]triphenylene is {final_answer} kcal/mol.")


solve_inversion_barrier()