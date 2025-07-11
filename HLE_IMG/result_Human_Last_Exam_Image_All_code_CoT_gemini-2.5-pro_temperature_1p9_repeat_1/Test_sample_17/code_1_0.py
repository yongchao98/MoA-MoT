import sympy as sp

def calculate_theoretical_stress_at_crack_tip():
    """
    This function demonstrates the theoretical stress at the tip of a sharp wedge
    using concepts from Linear Elastic Fracture Mechanics (LEFM).
    """

    # --- Step 1: Explain the theory ---
    print("This problem addresses stress concentration at a sharp notch (Point A).")
    print("In Linear Elastic Fracture Mechanics (LEFM), such a feature is modeled as a crack.")
    print("\nThe stress (sigma) directly ahead of a sharp crack tip is given by the equation:")
    # We will print the equation using its components
    equation_part1 = "sigma = K_I / sqrt("
    equation_part2 = "2 * pi * r"
    equation_part3 = ")"
    print(f"  {equation_part1}{equation_part2}{equation_part3}")

    print("\nWhere:")
    print("  K_I: The Stress Intensity Factor (a finite constant for a given load and geometry).")
    print("  pi: The mathematical constant pi.")
    print("  r: The distance from the crack tip.")
    print("  In the equation, the numbers are 2 and pi (approx 3.14159...).")


    # --- Step 2: Use symbolic math to find the stress AT the tip ---
    print("\nTo find the theoretical stress *at* the tip, we must calculate the limit of this equation as r -> 0.")

    # Define the symbols for our symbolic calculation
    K_I = sp.Symbol('K_I', real=True, positive=True) # K_I is a finite positive constant
    r = sp.Symbol('r', real=True, positive=True) # r is the distance

    # Define the expression for the stress
    # sp.pi represents the symbolic mathematical constant pi
    stress_expression = K_I / sp.sqrt(2 * sp.pi * r)

    # Calculate the limit as r approaches 0 from the positive side
    stress_at_tip = sp.limit(stress_expression, r, 0, dir='+')


    # --- Step 3: Print the conclusion ---
    print(f"\nUsing symbolic math, the limit of the expression '{stress_expression}' as r -> 0 is: {stress_at_tip}")
    print("\n----------------------------------------------------------------------------------------------------")
    print("Conclusion: The theoretical stress at the tip of an ideally sharp wedge is INFINITE.")
    print("Note: In real materials, this infinite stress is not physically possible. The material yields or fractures,")
    print("redistributing the stress over a small plastic zone at the tip.")
    print("----------------------------------------------------------------------------------------------------")

if __name__ == '__main__':
    calculate_theoretical_stress_at_crack_tip()
