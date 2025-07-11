def solve_projectile_explosion():
    """
    Calculates the landing position of the second fragment of an exploded projectile.
    """
    # Given value: horizontal distance from the gun to the highest point of elevation.
    I = 500  # in meters

    # The problem is solved using the principle of the center of mass (CoM).
    # The explosion is an internal force, so the trajectory of the CoM of the system is unchanged.
    # The CoM will land where the original projectile would have landed.
    
    # For a symmetric trajectory, the total range R would have been 2 * I.
    # R = 2 * I
    
    # The projectile splits into two equal fragments. Let their landing positions be x1 and x2.
    # The landing position of their CoM is R_cm = (x1 + x2) / 2.
    # Since R_cm = R, we get (x1 + x2) / 2 = 2 * I, which simplifies to:
    # x1 + x2 = 4 * I
    
    # It is given that "One of the fragments fell near the gun."
    # We interpret this to mean it landed at the gun's position, so x1 = 0.
    x1 = 0
    
    # We can now find the landing position of the second fragment, x2.
    
    print("The relationship between the landing positions of the two fragments (x1, x2) is derived from the conservation of the center of mass trajectory.")
    print("The derived equation is: x1 + x2 = 4 * I")
    print("\n--- Calculation ---")
    print(f"Given distance to the highest point, I = {I} m.")
    print(f"Landing position of the first fragment, x1 = {x1} m.")
    
    print("\nTo find the landing position of the second fragment (x2), we substitute the known values into the equation:")
    # Here we output each number in the final equation as requested.
    print(f"{x1} + x2 = 4 * {I}")
    
    # Perform the final calculation
    x2 = 4 * I - x1
    
    print(f"x2 = {x2}")
    
    print("\nTherefore, the maximum distance from the gun you can safely be is the landing spot of the second fragment.")
    print(f"Maximum safe distance = {x2} m.")

solve_projectile_explosion()
<<<2000>>>