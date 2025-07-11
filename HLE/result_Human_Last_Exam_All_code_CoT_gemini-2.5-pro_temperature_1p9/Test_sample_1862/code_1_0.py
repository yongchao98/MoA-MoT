def calculate_lithiation_time():
    """
    Calculates the lithiation time for a quarter-sphere LFP particle based on a reference full-sphere particle.
    """
    # Parameters for the first case (full sphere)
    d1_nm = 50  # Diameter in nm
    t1_s = 20   # Lithiation time in seconds

    # Parameters for the second case (quarter sphere)
    d2_nm = 250 # Diameter in nm

    # The characteristic diffusion length (L) is the radius of the particle.
    L1_nm = d1_nm / 2
    L2_nm = d2_nm / 2

    # The lithiation time (t) is proportional to the square of the diffusion length (L^2).
    # t2 / t1 = (L2 / L1)^2
    # t2 = t1 * (L2 / L1)^2
    t2_s = t1_s * (L2_nm / L1_nm)**2

    print("--- Calculation Steps ---")
    print(f"1. Initial particle (full sphere):")
    print(f"   - Diameter (d1): {d1_nm} nm")
    print(f"   - Radius (L1): {d1_nm} / 2 = {L1_nm} nm")
    print(f"   - Lithiation time (t1): {t1_s} s")
    print("\n2. Second particle (quarter sphere):")
    print(f"   - Diameter (d2): {d2_nm} nm")
    print(f"   - Radius (L2): {d2_nm} / 2 = {L2_nm} nm")
    print("\n3. Time calculation based on the relationship t ∝ L^2:")
    print(f"   t2 = t1 * (L2 / L1)^2")
    print(f"   t2 = {t1_s} s * ({L2_nm} nm / {L1_nm} nm)^2")
    print(f"   t2 = {t1_s} * ({L2_nm/L1_nm})**2")
    print(f"   t2 = {t1_s} * { (L2_nm/L1_nm)**2 }")
    print(f"   t2 = {t2_s} s")

    print("\n--- Final Answer ---")
    print(f"The calculated lithiation time for the quarter-sphere particle is {t2_s} seconds.")

    # Final answer in the required format
    print(f"\n<<<{t2_s}>>>")

calculate_lithiation_time()