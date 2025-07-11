import numpy as np

def generate_xrd_index_table():
    """
    Identifies prominent peaks from the CuO XRD pattern, calculates their
    d-spacing, and assigns the corresponding Miller indices (hkl).
    """

    # Standard X-ray wavelength for Cu K-alpha radiation (in Angstroms)
    lambda_wl = 1.5406

    # Prominent observed peaks (2-theta in degrees) from the provided XRD plot
    # and their corresponding Miller indices (hkl) for monoclinic CuO (Tenorite).
    # Data is based on standard reference patterns (e.g., JCPDS card no. 45-0937).
    peak_data = [
        (32.5, "(1 1 0)"),
        (35.5, "(-1 1 1)"),
        (38.7, "(1 1 1)"),
        (48.8, "(-2 0 2)"),
        (53.5, "(2 0 2)"),
        (58.3, "(-1 1 3)")
    ]

    # Print the table header
    print("Indexing of Prominent XRD Peaks for CuO Nanoparticles")
    print("-" * 55)
    print(f"{'2θ (Observed, deg)':<22} {'d-spacing (Å)':<20} {'Miller Index (h k l)':<20}")
    print("-" * 55)

    # Calculate d-spacing for each peak and print the table row
    for two_theta_deg, hkl in peak_data:
        # Convert 2-theta from degrees to theta in radians
        theta_rad = np.deg2rad(two_theta_deg / 2)
        
        # Calculate d-spacing using Bragg's Law: nλ = 2dsin(θ) -> d = λ / (2sin(θ)) for n=1
        d_spacing = lambda_wl / (2 * np.sin(theta_rad))
        
        # Print the results in a formatted table
        print(f"{two_theta_deg:<22.1f} {d_spacing:<20.4f} {hkl:<20}")
    
    print("-" * 55)

# Run the function to generate and print the table
generate_xrd_index_table()