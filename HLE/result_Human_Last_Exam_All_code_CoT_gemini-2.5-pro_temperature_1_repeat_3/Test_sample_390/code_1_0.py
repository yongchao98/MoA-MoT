import numpy as np

def solve_shape_problem():
    """
    Demonstrates that for orthogonal vectors y_i, the set S forms a simplex
    by finding the equation of the hyperplane it lies on.
    """
    # Define a set of 3 linearly independent and orthogonal vectors in R^3
    y1 = np.array([2.0, 0.0, 0.0])
    y2 = np.array([0.0, -3.0, 0.0])
    y3 = np.array([0.0, 0.0, 4.0])
    
    vectors = [y1, y2, y3]
    n = len(vectors)
    
    print("Given vectors:")
    for i, y in enumerate(vectors):
        print(f"y{i+1} = {y}")
    print("-" * 30)

    # Calculate the coefficients for the simplex equation
    # The coefficient for x_i is 1 / ||y_i||^2
    coeffs = [1.0 / np.linalg.norm(y)**2 for y in vectors]
    
    print("The shape S is a simplex defined by the linear equation:")
    equation_parts = []
    for i in range(n):
        equation_parts.append(f"({coeffs[i]:.4f}) * x{i+1}")
    equation = " + ".join(equation_parts) + " = 1"
    print(equation)
    print("\nWhich corresponds to:")
    equation_parts_frac = []
    for i, y in enumerate(vectors):
        equation_parts_frac.append(f"x{i+1} / {np.linalg.norm(y)**2:.1f}")
    equation_frac = " + ".join(equation_parts_frac) + " = 1"
    print(equation_frac)
    print("-" * 30)
    
    # Verification with a random vector s
    print("Verification:")
    # Generate a random unit vector s in the span of {y_i}
    # Since y_i are orthogonal and span R^3, we can generate any random unit vector in R^3
    s = np.random.randn(n)
    s = s / np.linalg.norm(s)
    print(f"Generated a random unit vector s = {s}")
    
    # Calculate the point (x1, x2, x3) in S
    x = np.array([np.abs(np.dot(y, s))**2 for y in vectors])
    print(f"Corresponding point in S: x = {x}")
    
    # Verify if it satisfies the equation
    total_sum = np.dot(coeffs, x)
    
    print("\nChecking if the point satisfies the equation:")
    check_equation_parts = []
    for i in range(n):
        check_equation_parts.append(f"({coeffs[i]:.4f}) * {x[i]:.4f}")
    check_equation = " + ".join(check_equation_parts)
    print(f"{check_equation} = {total_sum:.4f}")
    
    # The result should be very close to 1
    assert np.isclose(total_sum, 1.0)
    print("\nThe sum is indeed 1, confirming the point lies on the simplex.")

solve_shape_problem()
