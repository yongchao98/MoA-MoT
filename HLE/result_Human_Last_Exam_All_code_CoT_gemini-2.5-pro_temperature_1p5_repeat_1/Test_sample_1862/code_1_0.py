import math

def calculate_lithiation_time():
    """
    Calculates the lithiation time for a quarter-sphere LFP particle based on
    the given data for a smaller spherical particle.
    """
    # Parameters for the first case (spherical particle)
    d1_nm = 50.0  # Diameter in nm
    t1_s = 20.0   # Time in seconds

    # Parameters for the second case (quarter-sphere particle)
    d2_nm = 250.0 # Diameter in nm

    # The characteristic diffusion length (L) is the radius (r) of the particle.
    r1_nm = d1_nm / 2.0
    r2_nm = d2_nm / 2.0

    # The lithiation time (t) is proportional to the square of the diffusion length (r).
    # t_new / t_initial = (r_new / r_initial)^2
    # Therefore, t_new = t_initial * (r_new / r_initial)^2
    t2_s = t1_s * (r2_nm / r1_nm)**2

    print("The lithiation time (t) is proportional to the square of the particle's radius (r).")
    print("The relationship is: t_new = t_initial * (r_new / r_initial)²")
    print("\nGiven parameters:")
    print(f"Initial radius (r1): {d1_nm} nm / 2 = {r1_nm} nm")
    print(f"Initial time (t1): {t1_s} s")
    print(f"New radius (r2): {d2_nm} nm / 2 = {r2_nm} nm")
    print("\nFinal Equation:")
    # Print the equation with all numbers as requested
    print(f"New time = {t1_s} * ({r2_nm} / {r1_nm})²")
    print(f"New time = {t1_s} * {math.pow(r2_nm / r1_nm, 2)}")
    print(f"\nThe calculated lithiation time for the quarter-sphere particle is: {t2_s} seconds.")

calculate_lithiation_time()
<<<500>>>