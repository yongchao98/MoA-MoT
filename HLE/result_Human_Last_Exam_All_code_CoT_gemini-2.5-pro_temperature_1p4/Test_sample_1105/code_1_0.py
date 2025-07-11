import math

def calculate_max_distance():
    """
    Calculates the maximum distance the spaceship can achieve from the asteroid.
    
    The calculation is based on the work-energy theorem, leading to a quadratic equation
    for the maximum distance, l_max.
    """
    # Constants and initial conditions (using example values)
    # Gravitational constant
    G = 6.67430e-11  # m^3 kg^-1 s^-2
    # Mass of asteroid A
    m = 6e20  # kg (a large asteroid)
    # Mass of spaceship B
    M = 50000.0  # kg
    # Initial distance
    l0 = 20000.0  # m
    # Initial speed
    v0 = 200.0  # m/s
    # Applied force
    F = 150000.0 # N

    print("Given parameters:")
    print(f"  Gravitational constant (G): {G} m^3 kg^-1 s^-2")
    print(f"  Asteroid mass (m): {m} kg")
    print(f"  Spaceship mass (M): {M} kg")
    print(f"  Initial distance (l0): {l0} m")
    print(f"  Initial speed (v0): {v0} m/s")
    print(f"  Applied force (F): {F} N")
    print("-" * 30)

    # For a finite l_max to exist, the initial gravitational force must be greater
    # than the applied force. Otherwise, the spaceship accelerates indefinitely.
    if F >= (G * m * M) / (l0**2):
        print("The applied force is strong enough to cause indefinite acceleration.")
        print("The maximum distance is infinite.")
        return

    # From the work-energy theorem, we get the equation:
    # F*l_max^2 - C*l_max + G*m*M = 0
    # where C = F*l0 + (G*m*M)/l0 - 0.5*M*v0^2

    # Calculate the terms
    GmM = G * m * M
    Fl0 = F * l0
    GmM_div_l0 = GmM / l0
    Ki = 0.5 * M * v0**2

    # Calculate C
    C = Fl0 + GmM_div_l0 - Ki
    
    # Calculate the discriminant of the quadratic equation: D = C^2 - 4*F*(G*m*M)
    discriminant = C**2 - 4 * F * GmM

    if discriminant < 0:
        print("The initial velocity is too high for the spaceship to have a turning point.")
        print("The maximum distance is infinite.")
        return
        
    # Solve the quadratic equation for l_max using the '+' root for the maximum distance
    # l_max = (C + sqrt(discriminant)) / (2*F)
    l_max = (C + math.sqrt(discriminant)) / (2 * F)

    # Print the equation with the calculated numbers
    print("The maximum distance l_max is the solution to the equation:")
    print("F * l_max^2 - C * l_max + (G * m * M) = 0")
    print("Substituting the numerical values:")
    print(f"{F:.2f} * l_max^2 - ({C:.2f}) * l_max + ({GmM:.2f}) = 0\n")

    print("Solving using the quadratic formula l_max = (C + sqrt(C^2 - 4*F*GmM)) / (2*F):")
    print(f"l_max = ({C:.2f} + sqrt({C:.2f}^2 - 4 * {F:.2f} * {GmM:.2f})) / (2 * {F:.2f})")
    print(f"l_max = ({C:.2f} + sqrt({discriminant:.2f})) / {2 * F:.2f}")
    print(f"l_max = ({C:.2f} + {math.sqrt(discriminant):.2f}) / {2 * F:.2f}")
    print(f"l_max = {(C + math.sqrt(discriminant)):.2f} / {2 * F:.2f}\n")
    
    print("Result:")
    print(f"The maximum distance l_max is: {l_max:.2f} meters, or {l_max/1000:.2f} kilometers.")
    
    # The final answer in the requested format
    return l_max

# Execute the calculation
final_distance = calculate_max_distance()
# print(f"\n<<<Result: {final_distance}>>>") # Suppress this as per user request to use print inside the function.
# The user might be expecting a single literal answer at the very end. Let's provide that.
# The final result is already printed by the function. Let me formulate the final output string.
final_equation = f"l_{\\text{{max}}} = \\frac{{({F} \\cdot {l_0} + \\frac{{{G} \\cdot {m} \\cdot {M}}}{{{l_0}}} - \\frac{{1}}{{2}} \\cdot {M} \\cdot {v_0}^2) + \\sqrt{{({F} \\cdot {l_0} + \\frac{{{G} \\cdot {m} \\cdot {M}}}{{{l_0}}} - \\frac{{1}}{{2}} \\cdot {M} \\cdot {v_0}^2)^2 - 4 \\cdot {F} \\cdot {G} \\cdot {m} \\cdot {M}}}}}{{2 \\cdot {F}}}"
# I am not supposed to use this final format. I should just return the number.
# "directly return the answer with the format <<<answer content>>> at the end of your response"
# So I should print the number from the function call.
print(f"<<<{final_distance}>>>")