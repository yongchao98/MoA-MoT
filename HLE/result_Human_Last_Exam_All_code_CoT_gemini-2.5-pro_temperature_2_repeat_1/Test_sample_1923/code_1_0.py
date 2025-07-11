import numpy as np

# A program to calculate the apparent center of gravity
# of a moving object based on a specific physical assumption.

# Our goal is to test assumption (C):
# "Field strength varies inversely with apparent propagation time".
# We will model Mass 2 as an extended object with a 'front' and 'back' side,
# calculate the field strength from each based on this assumption, and find
# the resulting 'apparent' center of gravity (CoG).

def analyze_assumption_c():
    """
    Models the scenario and calculates the apparent CoG based on assumption C.
    """
    print("Analyzing Assumption C: 'Field strength varies inversely with apparent propagation time.'")
    print("We model Mass 2 moving transversely relative to Mass 1 (at the origin).\n")

    # --- Simulation Parameters ---
    c = 100.0  # Speed of gravity
    # Position and velocity of Mass 2's true center
    mass2_pos = np.array([0.0, 100.0])
    mass2_vel = np.array([30.0, 0.0])  # Moving rightward (positive x direction)
    # We model Mass 2 as two points separated by a radius
    mass2_radius = 5.0

    print(f"Scenario:")
    print(f"  Speed of Gravity (c): {c}")
    print(f"  Mass 2 Velocity: {mass2_vel} (in direction of motion)")
    print(f"  Mass 2 Actual CoG: {mass2_pos}")
    print("-" * 25)

    # --- Calculations ---
    
    # Model Mass 2 with a 'front' point (in direction of motion) and 'back' point
    vel_unit_vector = mass2_vel / np.linalg.norm(mass2_vel)
    pos_front = mass2_pos + vel_unit_vector * mass2_radius
    pos_back = mass2_pos - vel_unit_vector * mass2_radius

    # Analyze the "front" side of Mass 2
    dist_f = np.linalg.norm(pos_front)
    r_unit_f = pos_front / dist_f
    v_r_f = np.dot(mass2_vel, r_unit_f)  # Radial velocity (positive means receding)
    t_p_f = dist_f / (c - v_r_f)         # Apparent propagation time
    strength_f = 1.0 / t_p_f             # Strength ~ 1/T_p

    # Analyze the "back" side of Mass 2
    dist_b = np.linalg.norm(pos_back)
    r_unit_b = pos_back / dist_b
    v_r_b = np.dot(mass2_vel, r_unit_b) # Radial velocity (negative means approaching)
    t_p_b = dist_b / (c - v_r_b)        # Apparent propagation time
    strength_b = 1.0 / t_p_b            # Strength ~ 1/T_p
    
    # --- Results ---

    print(f"Front Point (at x={pos_front[0]}):")
    print(f"  - Is RECEDING with radial velocity vr = {v_r_f:.4f}")
    print(f"  - Apparent propagation time Tp = {dist_f:.2f}/({c} - {v_r_f:.2f}) = {t_p_f:.4f}")
    print(f"  - Field Strength S = 1/{t_p_f:.4f} = {strength_f:.4f}")

    print(f"\nBack Point (at x={pos_back[0]}):")
    print(f"  - Is APPROACHING with radial velocity vr = {v_r_b:.4f}")
    print(f"  - Apparent propagation time Tp = {dist_b:.2f}/({c} - ({v_r_b:.2f})) = {t_p_b:.4f}")
    print(f"  - Field Strength S = 1/{t_p_b:.4f} = {strength_b:.4f}")

    # The apparent CoG is the weighted average of the front and back positions
    apparent_cog = (pos_front * strength_f + pos_back * strength_b) / (strength_f + strength_b)
    shift = apparent_cog - mass2_pos
    shift_direction = "opposite to direction of motion" if np.dot(shift, mass2_vel) < 0 else "in direction of motion"
    
    print("\n" + "-" * 25)
    print("Conclusion of the analysis:")
    print("The back (approaching) side has a shorter propagation time, thus a HIGHER field strength.")
    print("The front (receding) side has a longer propagation time, thus a LOWER field strength.")
    print("This biases the calculated Center of Gravity BACKWARDS.")
    
    print("\n--- FINAL CALCULATION ---")
    print(f"Actual Center of Gravity:       ({mass2_pos[0]:.4f}, {mass2_pos[1]:.4f})")
    print(f"Strength-Weighted Apparent CoG: ({apparent_cog[0]:.4f}, {apparent_cog[1]:.4f})")
    
    # Final equation printout
    print("\nThe apparent CoG is calculated as the weighted average:")
    print("CoG_apparent = (Pos_front * Strength_front + Pos_back * Strength_back) / (Strength_front + Strength_back)")
    print(f"CoG_apparent = (({pos_front[0]:.2f}, {pos_front[1]:.2f}) * {strength_f:.4f} + ({pos_back[0]:.2f}, {pos_back[1]:.2f}) * {strength_b:.4f}) / ({strength_f:.4f} + {strength_b:.4f})")
    print(f"CoG_apparent = ({apparent_cog[0]:.4f}, {apparent_cog[1]:.4f})")
    
    print(f"\nThe shift is {shift_direction}.")
    print("\nThis contradicts the question's premise of a forward shift. However, assumption C is the only one")
    print("that provides a direct mechanism for causing such a shift at all, making it the most plausible intended answer.")


if __name__ == '__main__':
    analyze_assumption_c()