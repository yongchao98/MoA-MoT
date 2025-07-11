import math

def calculate_fusion_atoms():
    """
    Calculates the minimum number of titanium-50 atoms for a hypothetical
    sustained fusion reaction based on the Lawson criterion.
    """
    # Step 1: Define the constants and given parameters.

    # The Lawson criterion for ignition (sustained reaction) is approximately n*tau >= 5e20 s/m^3.
    # This is a benchmark for D-T fusion, used here for the hypothetical Ti-Ti reaction.
    lawson_ignition_product = 5e20  # units: s * m^-3

    # The particle confinement time is given as 'optimistic' 1 second.
    confinement_time_s = 1.0  # units: s

    # The reaction chamber is a cube with a side length of 10 cm.
    side_length_cm = 10.0  # units: cm

    # Step 2: Convert units for consistency.
    # Convert side length from centimeters to meters.
    side_length_m = side_length_cm / 100.0  # units: m

    # Step 3: Perform the calculations.

    # Calculate the volume of the cubic reaction chamber.
    # Volume = side_length^3
    volume_m3 = side_length_m ** 3

    # Calculate the required particle density (n) using the Lawson criterion.
    # n = (n*tau) / tau
    particle_density_per_m3 = lawson_ignition_product / confinement_time_s

    # Calculate the total number of atoms required to achieve this density in the given volume.
    # Total Atoms = density * volume
    total_atoms = particle_density_per_m3 * volume_m3

    # Step 4: Print the results, showing the final equation.
    print("--- Calculation for Minimum Required Atoms ---")
    print(f"Assumed Lawson Criterion for Ignition (n*tau): {lawson_ignition_product:.1e} s/m^3")
    print(f"Given Confinement Time (tau): {confinement_time_s} s")
    print(f"Given Chamber Side Length: {side_length_cm} cm ({side_length_m} m)")
    print("-" * 44)
    print(f"1. Chamber Volume = ({side_length_m} m)^3 = {volume_m3:.3f} m^3")
    print(f"2. Required Particle Density (n) = {lawson_ignition_product:.1e} / {confinement_time_s} = {particle_density_per_m3:.1e} particles/m^3")
    print("-" * 44)
    print("Final Equation: Total Atoms = Particle Density * Volume")
    # The final print statement shows each number in the final equation as requested.
    print(f"Total Atoms = {particle_density_per_m3:.1e} particles/m^3 * {volume_m3:.3f} m^3 = {total_atoms:.1e} atoms")

calculate_fusion_atoms()
print("\n<<<5e+17>>>")