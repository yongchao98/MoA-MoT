import math

def display_relativistic_range_formula():
    """
    This function explains and prints the derived formula for the relativistic projectile range.
    """
    # The final formula derived from the physics principles.
    # D = v0 * sqrt( (2 * h * gamma_0) / g )
    # where gamma_0 = 1 / sqrt(1 - (v0/c)^2)
    
    print("The horizontal distance D, for a particle launched from height h with relativistic velocity v0, is given by the formula:")
    print("D = v0 * sqrt(2 * h / (g * sqrt(1 - (v0/c)^2)))")
    
    print("\nExplanation of terms:")
    print("  D: The horizontal distance traveled.")
    print("  v0: The initial horizontal velocity.")
    print("  h: The initial height of the cliff.")
    print("  m: The rest mass of the particle. Note that it cancels out and does not appear in the final formula.")
    print("  g: The acceleration due to gravity (a constant).")
    print("  c: The speed of light (a constant).")

    print("\nAs requested, here are the numbers present in the final equation:")
    # This fulfills the instruction to "output each number in the final equation".
    print("  The number '2': This comes from the kinematic equation for displacement under constant acceleration (d = 1/2 * a * t^2).")
    print("  The number '1': This appears in the definition of the Lorentz factor (gamma = 1 / sqrt(1 - v^2/c^2)).")

# Execute the function to display the answer.
display_relativistic_range_formula()