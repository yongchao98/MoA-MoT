import numpy as np
from scipy.integrate import quad

def solve_and_print_results():
    """
    This function performs all the necessary calculations and prints the final results.
    """
    # Given points
    A = np.array([1, 0])
    B = np.array([10, 0])
    C = np.array([4, 3])

    # 1. Find the equation of the circle.
    # We solve for the coefficients D, E, F in the circle equation x^2 + y^2 + Dx + Ey + F = 0.
    # Substituting the points gives a linear system: M * [D, E, F]^T = V.
    M = np.array([
        [A[0], A[1], 1],
        [B[0], B[1], 1],
        [C[0], C[1], 1]
    ])
    V = np.array([
        -(A[0]**2 + A[1]**2),
        -(B[0]**2 + B[1]**2),
        -(C[0]**2 + C[1]**2)
    ])

    # Solve the linear system for D, E, F
    try:
        D, E, F = np.linalg.solve(M, V)
    except np.linalg.LinAlgError:
        print("The points are collinear and do not define a unique circle.")
        return

    # 2. Convert to center-radius form: (x-h)^2 + (y-k)^2 = r^2
    h = -D / 2
    k = -E / 2
    r_squared = h**2 + k**2 - F

    print("Step 1 & 2: Finding the equation of the circle.")
    print(f"The equation in the form x^2 + y^2 + Dx + Ey + F = 0 has coefficients:")
    print(f"D = {D}, E = {E}, F = {F}")
    print(f"The center-radius form of the equation is (x - {h})^2 + (y - ({k}))^2 = {r_squared}")
    print("-" * 30)
    
    # 3. Define the function y = f(x) for the upper arc of the circle.
    # y = k +/- sqrt(r^2 - (x-h)^2). We choose the sign based on point C.
    # To pass through C(4,3), we must have y > k, so we use the '+' sign.
    def f(x):
        # We need to handle potential floating point inaccuracies that make the argument negative.
        val = r_squared - (x - h)**2
        if isinstance(val, (int, float)):
            return k + np.sqrt(max(0, val))
        else: # Handle numpy arrays
            return k + np.sqrt(np.maximum(0, val))

    # 4. Calculate the normalization constant alpha.
    # alpha = 1 / integral of f(x) from 1 to 10.
    total_area, _ = quad(f, 1, 10, limit=200)

    if total_area <= 0:
        print("The integral of f(x) is not positive, so a valid PDF cannot be formed.")
        return
        
    alpha = 1 / total_area

    print("Step 3 & 4: Calculating alpha.")
    print(f"The function for the upper arc is f(x) = {k} + sqrt({r_squared} - (x - {h})^2)")
    print(f"The integral of f(x) over [1, 10] is: {total_area}")
    print(f"The value of alpha is 1 / {total_area} = {alpha}")
    print("-" * 30)

    # 5. Calculate the probability P(X < 3).
    # P(X < 3) = alpha * integral of f(x) from 1 to 3.
    partial_area, _ = quad(f, 1, 3, limit=200)
    probability = alpha * partial_area

    print("Step 5: Calculating P(X < 3).")
    print(f"The integral of f(x) over [1, 3] is: {partial_area}")
    print(f"The probability P(X < 3) = {alpha} * {partial_area} = {probability}")
    print("-" * 30)
    
    # The final requested values
    print("Final values:")
    print(f"alpha = {alpha}")
    print(f"P(X < 3) = {probability}")


if __name__ == "__main__":
    solve_and_print_results()