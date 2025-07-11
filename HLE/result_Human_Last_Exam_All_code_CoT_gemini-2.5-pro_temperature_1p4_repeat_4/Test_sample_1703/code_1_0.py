import math

def calculate_area_coefficients():
    """
    This function calculates the coefficients for the area formula of the triangle T(t).
    The area is a quadratic function of time, Area(t) = a*t^2 + b.
    """
    
    # The coefficient 'a' for the t^2 term is derived from the geometry and kinematics.
    # a = (3 * sqrt(3)) / 4
    a = (3 * math.sqrt(3)) / 4
    
    # The constant term 'b' is the area at t=0.
    # b = (225 * sqrt(3)) / 4
    b = (225 * math.sqrt(3)) / 4
    
    return a, b

def main():
    """
    Calculates and prints the final equation for the area of triangle T(t).
    """
    
    a_coeff, b_coeff = calculate_area_coefficients()
    
    print("The area of the triangle T(t) is a quadratic function of time t.")
    print("The formula is in the form: Area(t) = a * t^2 + b\n")
    
    print("The numbers (coefficients) for this equation are:")
    print(f"a = {a_coeff}")
    print(f"b = {b_coeff}\n")
    
    print("The final equation with the calculated numbers is:")
    print(f"Area(t) = {a_coeff} * t^2 + {b_coeff}\n")
    
    print("For precision, the symbolic form of the equation is:")
    print("Area(t) = (3*sqrt(3)/4) * t^2 + (225*sqrt(3)/4)")
    print("\nThis simplifies to:")
    print("Area(t) = (3*sqrt(3)/4) * (t^2 + 75)")

if __name__ == "__main__":
    main()