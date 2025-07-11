def calculate_lithiation_time():
    """
    Calculates the lithiation time for a quarter-sphere LFP particle based on
    a surface reaction rate-limited model.
    """
    # --- Parameters for the initial spherical particle ---
    d1_nm = 50.0
    r1_nm = d1_nm / 2.0
    t1_s = 20.0

    # --- Parameters for the final quarter-sphere particle ---
    d2_nm = 250.0
    r2_nm = d2_nm / 2.0

    # --- Calculation ---
    # The lithiation time (t) is proportional to the particle's volume-to-surface-area ratio.
    # The relationship is: t2 = t1 * (r2 / r1) * 0.5
    t2_s = t1_s * (r2_nm / r1_nm) * 0.5

    # --- Output ---
    print("Based on a surface reaction-limited model, the final lithiation time (t2) is calculated.")
    print("\nThe governing equation is:")
    print("t2 = t1 * (radius_2 / radius_1) * 0.5")
    print("\nWhere:")
    print(f"t1 (initial time) = {t1_s} s")
    print(f"radius_1 (initial radius) = {r1_nm} nm")
    print(f"radius_2 (final radius) = {r2_nm} nm")
    
    print("\nFinal Calculation:")
    # Printing each number in the final equation
    print(f"t2 = {t1_s} * ({r2_nm} / {r1_nm}) * 0.5")
    print(f"t2 = {t1_s} * {r2_nm/r1_nm} * 0.5")
    print(f"t2 = {t2_s} s")

# Execute the function to print the result
calculate_lithiation_time()
<<<50.0>>>