import math

def analyze_proton_beam(total_laser_energy=10.0, gaussian_radius=1.0, oam_outer_radius=2.0, oam_inner_radius=1.5):
    """
    This function conceptually models the effect of a laser beam's shape on proton acceleration.
    It compares a standard Gaussian beam with a beam carrying Orbital Angular Momentum (OAM).
    """

    # --- Step 1: Calculate the area and intensity for a standard Gaussian beam ---
    # A Gaussian beam focuses energy to a small spot.
    # Area = pi * r^2
    area_gaussian = math.pi * gaussian_radius**2
    # Intensity = Energy / Area
    intensity_gaussian = total_laser_energy / area_gaussian

    # --- Step 2: Calculate the area and intensity for an OAM beam ---
    # An OAM beam has a doughnut-shaped (annular) profile, spreading the energy out.
    # Area = pi * (r_outer^2 - r_inner^2)
    area_oam = math.pi * (oam_outer_radius**2 - oam_inner_radius**2)
    intensity_oam = total_laser_energy / area_oam

    # --- Step 3: Relate laser intensity to proton energy ---
    # A common scaling law states that maximum proton energy is proportional
    # to the square root of the peak laser intensity: E_proton ∝ sqrt(I)
    # We'll use a proportionality constant of 1 for this relative comparison.
    k = 1.0
    energy_proton_gaussian_eq = f"{k} * sqrt({intensity_gaussian:.2f})"
    energy_proton_oam_eq = f"{k} * sqrt({intensity_oam:.2f})"
    
    energy_proton_gaussian = k * math.sqrt(intensity_gaussian)
    energy_proton_oam = k * math.sqrt(intensity_oam)


    # --- Step 4: Print the analysis and conclusion ---
    print("Analysis of Laser to Proton Conversion")
    print("=" * 40)
    
    # Conclusion on Energy
    print("\n[Proton Energy Analysis]")
    print("A key principle is that higher laser intensity leads to higher proton energy.")
    print(f"  - Gaussian Beam Intensity: {total_laser_energy} / {area_gaussian:.2f} = {intensity_gaussian:.2f} units")
    print(f"  - OAM Beam Intensity:      {total_laser_energy} / {area_oam:.2f} = {intensity_oam:.2f} units")
    print("\nThe OAM beam's larger area results in lower intensity.")
    print("Using the scaling law 'E_proton ∝ sqrt(Intensity)':")
    print(f"  - Relative Gaussian Proton Energy = sqrt({intensity_gaussian:.2f}) = {energy_proton_gaussian:.2f} units")
    print(f"  - Relative OAM Proton Energy      = sqrt({intensity_oam:.2f}) = {energy_proton_oam:.2f} units")
    print("Conclusion 1: The proton energy DECREASES.")

    # Conclusion on Dispersion
    print("\n[Proton Beam Collimation Analysis]")
    print("The OAM beam's ring-shaped profile creates a ring-shaped accelerating field.")
    print("This field exerts an outward force on the protons as they are accelerated forward.")
    print("Conclusion 2: This leads to the DISPERSION of the proton beam.")

    print("\n" + "=" * 40)
    print("Final Result: The proton beam experiences DISPERSION and PROTON ENERGY DECREASES.")
    print("=" * 40)


if __name__ == '__main__':
    analyze_proton_beam()
<<<C>>>