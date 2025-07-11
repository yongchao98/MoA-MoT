import math

def calculate_lithiation_time():
    """
    Calculates the lithiation time for a quarter-sphere LFP particle based on
    the known lithiation time for a smaller, full-sphere particle.
    """
    # Given parameters for the first case (full sphere)
    d1_nm = 50.0  # Diameter in nm
    t1_s = 20.0   # Lithiation time in seconds

    # Given parameters for the second case (quarter sphere)
    d2_nm = 250.0 # Diameter of the original sphere in nm

    # --- Calculations ---

    # Step 1: Calculate the radii
    r1_nm = d1_nm / 2.0
    r2_nm = d2_nm / 2.0

    # Step 2: Use the derived relationship to find the new lithiation time (t2).
    # The relationship is: t2 = t1 * (2/3) * (r2 / r1)^2
    
    # Calculate the ratio of the radii and its square
    radius_ratio = r2_nm / r1_nm
    radius_ratio_squared = radius_ratio ** 2
    
    # Calculate the final time
    t2_s = t1_s * (2.0 / 3.0) * radius_ratio_squared

    # --- Output Results ---
    
    print("The physical model for diffusion-limited full lithiation gives the relationship:")
    print("t_new = t_initial * (2/3) * (r_new / r_initial)^2\n")
    print("Initial parameters:")
    print(f"Particle 1 (full sphere): Diameter = {d1_nm} nm, Radius (r1) = {r1_nm} nm, Time (t1) = {t1_s} s")
    print(f"Particle 2 (quarter sphere): Original Diameter = {d2_nm} nm, Radius (r2) = {r2_nm} nm\n")
    print("Calculation:")
    print(f"Lithiation Time = {t1_s} s * (2/3) * ({r2_nm} nm / {r1_nm} nm)^2")
    print(f"Lithiation Time = {t1_s} s * (2/3) * ({radius_ratio})^2")
    print(f"Lithiation Time = {t1_s} s * (2/3) * {radius_ratio_squared}")
    print(f"Lithiation Time = {(t1_s * 2 * radius_ratio_squared):.1f} / 3 s")
    print(f"Final Lithiation Time = {t2_s:.2f} s")

calculate_lithiation_time()
<<<333.33>>>