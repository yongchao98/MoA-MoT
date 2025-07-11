import math

def calculate_lithiation_time():
    """
    Calculates the lithiation time for a quarter-sphere LFP particle based on
    the lithiation time of a smaller spherical particle.
    """
    # Initial case parameters
    d1_nm = 50  # Initial particle diameter in nm
    t1_s = 20   # Initial lithiation time in seconds

    # Final case parameters
    d2_nm = 250 # Final particle diameter in nm

    # Calculate radii from diameters
    r1_nm = d1_nm / 2
    r2_nm = d2_nm / 2

    # The lithiation time (t) is governed by solid-state diffusion and is
    # proportional to the square of the diffusion length (the particle radius, r).
    # Relationship: t_final / t_initial = (r_final / r_initial)^2
    # This means: t_final = t_initial * (r_final / r_initial)^2

    t2_s = t1_s * (r2_nm / r1_nm)**2

    print("Step 1: Identify the relationship between lithiation time and particle size.")
    print("Lithiation time (t) is proportional to the square of the particle radius (r): t ∝ r².")
    print("\nStep 2: Set up the calculation based on the given values.")
    print(f"Initial Time (t_initial): {t1_s} s")
    print(f"Initial Radius (r_initial): {d1_nm} nm / 2 = {r1_nm} nm")
    print(f"Final Radius (r_final): {d2_nm} nm / 2 = {r2_nm} nm")
    
    print("\nStep 3: Calculate the final lithiation time (t_final) using the formula.")
    print("Formula: t_final = t_initial * (r_final / r_initial)²")
    
    # Print the final equation with the numbers plugged in
    print("\nCalculation with values:")
    print(f"t_final = {t1_s} * ({r2_nm} / {r1_nm})²")
    ratio = r2_nm / r1_nm
    print(f"t_final = {t1_s} * ({ratio})²")
    print(f"t_final = {t1_s} * {ratio**2}")
    
    print(f"\nThe new lithiation time is: {int(t2_s)} seconds.")
    print("\nNote: The change in shape from a full sphere to a quarter sphere does not affect the characteristic diffusion length (the radius), which is the limiting factor for complete lithiation time in this idealized scenario.")

calculate_lithiation_time()
<<<500>>>