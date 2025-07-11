import math

def calculate_max_distance():
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid
    based on a set of given physical parameters.
    """
    # Problem parameters (chosen to illustrate a non-escape scenario)
    m = 1.0e18      # Mass of the asteroid (kg)
    M = 1.0e5       # Mass of the spaceship (kg)
    l0 = 10000.0    # Initial distance (m)
    v0 = 50.0       # Initial velocity away from the asteroid (m/s)
    F = 1000.0      # Applied constant force on the spaceship (N)
    G = 6.67430e-11 # Gravitational constant in N m^2/kg^2

    print("--- Input Parameters ---")
    print(f"Asteroid mass (m): {m:.4g} kg")
    print(f"Spaceship mass (M): {M:.4g} kg")
    print(f"Initial distance (l0): {l0:.4g} m")
    print(f"Initial speed (v0): {v0:.4g} m/s")
    print(f"Applied force (F): {F:.4g} N\n")

    # The problem can be solved by finding the turning points (where velocity is zero).
    # This leads to a quadratic equation derived from the work-energy theorem:
    # a*l^2 + b*l + c = 0, where the coefficients are determined below.
    
    # --- Step 1: Calculate the coefficients of the quadratic equation ---
    # The full equation is: F*l^2 - (F*l0 + G*M*m/l0 - 0.5*M*v0^2)*l + G*M*m = 0
    
    a = F
    GmM = G * M * m
    C = F * l0 + GmM / l0 - 0.5 * M * v0**2
    b = -C
    c = GmM

    print("--- Step 2: Set up the quadratic equation for turning points (l) ---")
    print("Equation form: a*l^2 + b*l + c = 0")
    print(f"a = F = {a:.4e}")
    term_F_l0 = F * l0
    term_GmM_l0 = GmM / l0
    term_KE = 0.5 * M * v0**2
    print(f"b = -(F*l0 + G*M*m/l0 - 0.5*M*v0^2) = -({term_F_l0:.4e} + {term_GmM_l0:.4e} - {term_KE:.4e}) = {b:.4e}")
    print(f"c = G*M*m = {c:.4e}")
    print(f"Final Equation: ({a:.4e})*l^2 + ({b:.4e})*l + ({c:.4e}) = 0\n")

    # --- Step 3: Solve the quadratic equation ---
    discriminant = b**2 - 4 * a * c

    if discriminant < 0:
        print("Result: The discriminant is negative. The spaceship escapes to infinity.")
        final_answer = float('inf')
    else:
        sqrt_discriminant = math.sqrt(discriminant)
        
        # The two roots (turning points) are l1 and l2
        # We ensure l1 <= l2 by convention
        l1 = (-b - sqrt_discriminant) / (2 * a)
        l2 = (-b + sqrt_discriminant) / (2 * a)

        print("--- Step 4: Find the turning points using the quadratic formula ---")
        print(f"Formula: l = (-b ± sqrt(b^2 - 4ac)) / 2a")
        print(f"Discriminant D = b^2 - 4ac = ({b:.4e})^2 - 4*({a:.4e})*({c:.4e}) = {discriminant:.4e}")
        print(f"The turning points are:")
        print(f"l1 = ({-b:.4e} - {sqrt_discriminant:.4e}) / (2 * {a:.4e}) = {l1:.5g} m")
        print(f"l2 = ({-b:.4e} + {sqrt_discriminant:.4e}) / (2 * {a:.4e}) = {l2:.5g} m\n")
        
        # --- Step 5: Determine the correct l_max ---
        print(f"--- Step 5: Determine the maximum distance from the turning points ---")
        print(f"Initial distance l0 = {l0:.5g} m.")

        if l0 > l2:
            print("l0 is greater than the outer turning point. The spaceship has enough energy to escape.")
            final_answer = float('inf')
        elif l0 < l1:
            print("l0 is less than the inner turning point. The spaceship moves away and stops at the first turning point it encounters.")
            final_answer = l1
            print(f"Result: l_max = l1 = {final_answer:.5g} m.")
        else: # l1 <= l0 <= l2
            print("l0 is between the turning points. The spaceship moves away towards the outer turning point.")
            final_answer = l2
            print(f"Result: l_max = l2 = {final_answer:.5g} m.")

    # Print the final answer in the required format
    if isinstance(final_answer, float):
        print(f"\n<<<{final_answer}>>>")

if __name__ == "__main__":
    calculate_max_distance()