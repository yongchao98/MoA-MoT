import math

def calculate_max_distance():
    """
    Calculates the maximum distance a spaceship can travel from an asteroid
    under gravity and an additional constant force.
    """
    # --- Input Parameters (using example values for demonstration) ---
    # Gravitational constant in m^3 kg^-1 s^-2
    G = 6.67430e-11
    # Mass of the asteroid in kg (e.g., a large asteroid like Vesta)
    m = 5.25e20
    # Mass of the spaceship in kg
    M = 5.0e4
    # Initial distance in meters (e.g., 1000 km)
    l_0 = 1.0e6
    # Initial speed in m/s (e.g., 200 m/s)
    v_0 = 200.0
    # Additional applied force from thrusters in Newtons
    F = 40.0

    print("This script calculates the maximum distance l_max based on the work-energy theorem.")
    print("The theorem leads to the following quadratic equation for l_max:")
    print("a * l_max^2 + b * l_max + c = 0\n")
    print("Where:")
    print("a = F")
    print("b = -(F*l_0 + G*M*m/l_0 - 0.5*M*v_0^2)")
    print("c = G*M*m\n")

    # --- Calculation Steps ---
    
    # 1. Calculate the coefficients of the quadratic equation
    a = F
    
    # This term represents the total initial energy configuration translated into the quadratic coefficient.
    energy_term = F * l_0 + (G * M * m) / l_0 - 0.5 * M * v_0**2
    b = -energy_term
    
    c = G * M * m

    # 2. Calculate the discriminant (the part under the square root in the quadratic formula)
    discriminant = b**2 - 4 * a * c

    # 3. Check if a real solution exists
    if discriminant < 0:
        print("Result: The initial speed and applied force are too high.")
        print("The spaceship will escape the asteroid's gravity and never stop.")
        print("There is no finite maximum distance (l_max -> infinity).")
        return

    # 4. Solve for the roots using the quadratic formula
    # The larger root corresponds to the maximum distance.
    l_max = (-b + math.sqrt(discriminant)) / (2 * a)

    # --- Output the Final Equation and Result ---
    
    print("The solution for l_max is the larger root of the quadratic formula:")
    print("l_max = (-b + sqrt(b^2 - 4ac)) / (2a)\n")
    
    print("Plugging in the numbers for the coefficients:")
    print(f"a = {a}")
    print(f"b = -({F}*{l_0} + ({G}*{M}*{m})/{l_0} - 0.5*{M}*{v_0**2}) = {b:.4e}")
    print(f"c = {G}*{M}*{m} = {c:.4e}\n")

    print("The final equation for l_max with the calculated numbers is:")
    print(f"l_max = ( -({b:.4e}) + sqrt( ({b:.4e})^2 - 4*({a})*({c:.4e}) ) ) / ( 2*({a}) )")
    print(f"l_max = ( {(-b):.4e} + sqrt( {discriminant:.4e} ) ) / ( {2*a} )")
    print(f"l_max = ( {(-b):.4e} + {math.sqrt(discriminant):.4e} ) / ( {2*a} )")

    # Print the final result
    print(f"\n<<<l_max = {l_max:.4e} meters>>>")

if __name__ == '__main__':
    calculate_max_distance()