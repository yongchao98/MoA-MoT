import sympy
from sympy import pi, integrate, Symbol, pretty_print, Eq

def calculate_enclosed_volume():
    """
    This function calculates the volume enclosed by a cone and an ellipsoid
    by following the method of slicing.
    """
    # Define the variable of integration
    y = Symbol('y')

    # Define the square of the radii for the cone and the ellipsoid as functions of y
    Rc_sq = ((4 - y) / 2)**2
    Re_sq = 3 * (1 - y**2 / 4)

    print("Step 1: Calculate the volume of the top part of the region (2 <= y <= 4).")
    print("This volume is bounded only by the cone.")
    # Volume V1 is the integral from y=2 to y=4 of the cone's cross-sectional area
    integrand1 = pi * Rc_sq
    print("The integral for the volume V1 is:")
    pretty_print(Eq(Symbol('V1'), integrate(integrand1, (y, 2, 4))))
    
    # Calculate the value of V1
    V1 = integrate(integrand1, (y, 2, 4))
    print(f"\nThe value of the volume of the top part is: V1 = {V1}\n")

    print("-" * 30)

    print("Step 2: Calculate the volume of the bottom part of the region (1 <= y <= 2).")
    print("This volume is the region between the cone and the ellipsoid.")
    # The integrand for V2 is the area of the annulus
    integrand2 = pi * (Rc_sq - Re_sq)
    print("The integral for the volume V2 is:")
    pretty_print(Eq(Symbol('V2'), integrate(integrand2, y)))
    
    # Let's simplify the integrand to make the integration clearer
    simplified_integrand2 = sympy.simplify(integrand2)
    print("\nSimplified integrand for V2:")
    pretty_print(simplified_integrand2)
    
    print("\nCalculating the definite integral for V2 from y=1 to y=2:")
    V2 = integrate(simplified_integrand2, (y, 1, 2))
    pretty_print(Eq(Symbol('V2'), integrate(simplified_integrand2, (y, 1, 2))))
    print(f"\nThe value of the volume of the bottom part is: V2 = {V2}\n")

    print("-" * 30)
    
    # Total volume is the sum of V1 and V2
    V_total = V1 + V2
    print("Step 3: Calculate the total volume.")
    print("The total enclosed volume is the sum of V1 and V2.")
    print(f"Total Volume V = V1 + V2")
    # Using the pretty_print to show the equation with symbols
    equation_str = f"{V1} + {V2} = {V_total}"
    print(f"V = {equation_str}")
    

if __name__ == '__main__':
    calculate_enclosed_volume()