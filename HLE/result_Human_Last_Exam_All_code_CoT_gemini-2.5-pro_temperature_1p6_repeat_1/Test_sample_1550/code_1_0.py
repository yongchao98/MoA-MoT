import sympy

def find_intersection_points():
    """
    Solves a system of two quadratic equations to find the number of real intersection points.
    """
    # 1. Define x and y as real symbolic variables
    x, y = sympy.symbols('x y', real=True)

    # 2. Define the two quadric equations from the problem statement
    # Equation 1: 164*x^2 - 216*x*y + 72*y^2 - 16*x + 31 = 0
    eq1 = 164*x**2 - 216*x*y + 72*y**2 - 16*x + 31

    # Equation 2: 864*x^2 - 1056*x*y + 324*y^2 - 560*x + 324*y + 149 = 0
    eq2 = 864*x**2 - 1056*x*y + 324*y**2 - 560*x + 324*y + 149

    # 3. Solve the system of equations for real (x, y) pairs
    solutions = sympy.nonlinsolve([eq1, eq2], [x, y])

    # 4. Count the number of solutions found
    # Since we declared x and y as real, nonlinsolve will only return real solutions.
    # The number of solutions is the length of the resulting FiniteSet.
    count = len(solutions)

    print(f"The system of equations has {count} real intersection point(s).")
    
    # For each solution found, print its coordinates
    if count > 0:
        print("The intersection points are:")
        for i, sol in enumerate(solutions):
            # Evaluate to floating point for clean printing if necessary, but Fraction is more exact.
            print(f"  Point {i+1}: (x={sol[0]}, y={sol[1]})")

    # 5. Fulfill the "output each number in the final equation" requirement
    # Let C be the count. We form the equation 1*x - C = 0.
    # The numbers in this equation are the coefficients 1 and -C.
    print("\nRepresenting the count C as the equation: 1*x - C = 0")
    print("The numbers in this final equation are:")
    
    # The coefficients are 1 and -count
    coeff1 = 1
    coeff2 = -count
    
    print(coeff1)
    print(coeff2)


if __name__ == '__main__':
    find_intersection_points()