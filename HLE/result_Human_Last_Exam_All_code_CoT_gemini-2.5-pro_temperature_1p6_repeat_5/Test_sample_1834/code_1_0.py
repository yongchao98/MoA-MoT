import numpy as np

def calculate_magnetic_field():
    """
    Calculates the magnitude of the magnetic field at a point P(1,-1,0)
    due to two infinite perpendicular wires carrying current I.
    """
    # --- Constants and Parameters ---
    # The permeability of free space (mu_0) in Tesla*meter/Ampere
    mu_0 = 4 * np.pi * 1e-7

    # We assume a current I = 1.0 A for a numerical result.
    # The final magnetic field magnitude is directly proportional to I.
    I = 1.0

    # The point P where the magnetic field is to be calculated
    P = np.array([1.0, -1.0, 0.0])

    print("--- Calculating Magnetic Field ---")
    print(f"Point of interest P: {P}")
    print(f"Current in each wire I: {I} A\n")

    # --- Contribution from Wire 1 (on x-axis) ---
    # Perpendicular distance r1 from the x-axis to point P
    r1 = np.sqrt(P[1]**2 + P[2]**2)
    # Magnitude of the magnetic field from Wire 1
    B1_mag = (mu_0 * I) / (2 * np.pi * r1)
    # By the right-hand rule, the field direction is -z.

    print("--- Wire 1 (on x-axis) ---")
    print(f"Perpendicular distance r1 = {r1:.2f} m")
    print(f"Magnitude B1 = {B1_mag:.2e} T (in the -z direction)")

    # --- Contribution from Wire 2 (on y-axis) ---
    # Perpendicular distance r2 from the y-axis to point P
    r2 = np.sqrt(P[0]**2 + P[2]**2)
    # Magnitude of the magnetic field from Wire 2
    B2_mag = (mu_0 * I) / (2 * np.pi * r2)
    # By the right-hand rule, the field direction is also -z.
    
    print("\n--- Wire 2 (on y-axis) ---")
    print(f"Perpendicular distance r2 = {r2:.2f} m")
    print(f"Magnitude B2 = {B2_mag:.2e} T (in the -z direction)")

    # --- Total Magnetic Field ---
    # Since both fields point in the same direction, the total magnitude
    # is the sum of the individual magnitudes.
    B_total_mag = B1_mag + B2_mag

    print("\n--- Total Magnetic Field ---")
    print("The total field is the sum of the fields from each wire.")
    print("Since both fields point in the same direction, we sum their magnitudes.")
    # The final equation with each number printed out
    print(f"|B_total| = |B1| + |B2| = {B1_mag:.2e} + {B2_mag:.2e} = {B_total_mag:.2e} T")
    print("\nFinal numerical answer for the magnitude (assuming I=1A):")
    print(B_total_mag)

if __name__ == '__main__':
    calculate_magnetic_field()
