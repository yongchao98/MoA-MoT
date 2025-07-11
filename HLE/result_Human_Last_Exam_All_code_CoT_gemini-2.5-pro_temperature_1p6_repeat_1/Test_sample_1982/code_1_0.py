import sys

def solve_control_problem():
    """
    Solves for the control variable u1 based on the given matrix equation and parameters.
    """
    # The matrix equation is:
    # [[0, 1], [0, 0]] * [[x11, x12], [x21, x22]] + [[x11, x12], [x21, x22]] * [[1, 0], [0, 1]] = 
    # = [[1, 0], [0, 1]] + [[c1, 0], [0, c2]] * [[u1, 0], [0, u2]]

    # Expanding the left side:
    # [[x21, x22], [0, 0]] + [[x11, x12], [x21, x22]] = [[x11 + x21, x12 + x22], [x21, x22]]

    # Expanding the right side:
    # [[1, 0], [0, 1]] + [[c1*u1, 0], [0, c2*u2]] = [[1 + c1*u1, 0], [0, 1 + c2*u2]]

    # Equating the two sides gives four scalar equations:
    # 1) x11 + x21 = 1 + c1 * u1
    # 2) x12 + x22 = 0
    # 3) x21 = 0
    # 4) x22 = 1

    # From equation (3), we get x21 = 0.
    # Substituting x21 = 0 into equation (1):
    # x11 + 0 = 1 + c1 * u1
    # This simplifies to x11 = 1 + c1 * u1.

    # Solving for u1 gives:
    # u1 = (x11 - 1) / c1

    # The problem provides c1 = 10^4 and alpha1. We assume x11 = alpha1.
    print("Deriving the formula for u_1 from the matrix equation:")
    print("1. The (2,1) element of the equation gives: x_21 = 0.")
    print("2. The (1,1) element gives: x_11 + x_21 = 1 + c_1 * u_1.")
    print("3. Substituting x_21 = 0 gives: x_11 = 1 + c_1 * u_1.")
    print("4. Solving for u_1 yields the formula: u_1 = (x_11 - 1) / c_1.\n")
    
    print("Assuming x_11 corresponds to the given parameter α_1.\n")

    # Given parameters
    c1 = 10**4
    # Note: Python's integers handle arbitrary size, so we can calculate alpha1 directly.
    alpha1 = (1 + 10**5)**6 * (1 - 10**5 + 10**10)

    print("Given values:")
    print(f"c_1 = {c1}")
    print(f"α_1 = (1 + 10^5)^6 * (1 - 10^5 + 10^10)")
    print(f"Calculated value of α_1 = {alpha1}\n")
    
    # We can now substitute the values into our formula for u_1.
    print("Substituting the numerical values into the formula for u_1:")
    print(f"u_1 = ({alpha1} - 1) / {c1}\n")
    
    # Calculate u1. The result is an integer, so we use integer division //.
    if (alpha1 - 1) % c1 != 0:
        print("Warning: The result is not a perfect integer.", file=sys.stderr)
        
    u1 = (alpha1 - 1) // c1
    
    print("Final result:")
    print(f"u_1 = {u1}")

solve_control_problem()