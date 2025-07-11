import numpy as np

def generate_cuo_xrd_table():
    """
    Calculates and prints a table of indexed XRD peaks for monoclinic CuO (tenorite),
    correlating observed peaks with theoretical values.
    """
    # --- Standard Data for CuO (Tenorite, monoclinic) ---
    # Source: JCPDS Card No. 45-0937 or similar crystallographic databases.
    a = 4.6837  # Lattice parameter a in Angstroms
    b = 3.4226  # Lattice parameter b in Angstroms
    c = 5.1288  # Lattice parameter c in Angstroms
    beta_deg = 99.54   # Beta angle in degrees
    beta_rad = np.radians(beta_deg)

    # --- X-ray source ---
    # Standard assumption is Copper K-alpha radiation
    lambda_val = 1.5406  # Wavelength in Angstroms

    print("XRD Peak Indexing for Monoclinic CuO (Tenorite)")
    print(f"Calculation based on Cu Kα radiation (λ = {lambda_val} Å) and lattice parameters:")
    print(f"a = {a} Å, b = {b} Å, c = {c} Å, β = {beta_deg}°\n")

    # Header for the results table
    header = f"{'2θ (Approx. From Graph) [°]':<30} {'Miller Indices (h k l)':<25} {'d-spacing (Calculated) [Å]':<30} {'2θ (Calculated) [°]':<20}"
    print(header)
    print("-" * len(header))

    # List of prominent (hkl) planes for CuO identified from the graph and literature
    # Each tuple contains: (approximate 2-theta from graph, (h, k, l) Miller indices)
    peaks_data = [
        (32.5, (1, 1, 0)),
        (35.5, (1, 1, -1)),
        (38.7, (1, 1, 1)),
        (48.7, (2, 0, -2)),
        (53.5, (2, 0, 2)),
        (58.3, (1, 1, -3)),
        (61.5, (-3, 1, 1)),
    ]

    # Loop through each identified peak to calculate its theoretical position
    for approx_2theta, (h, k, l) in peaks_data:
        # Monoclinic d-spacing formula: 1/d² = (1/sin²β) * (h²/a² + k²sin²β/b² + l²/c² - 2hl*cosβ/(ac))
        sin_beta_sq = np.sin(beta_rad)**2
        cos_beta = np.cos(beta_rad)
        
        # Calculate the inverse of d-squared
        d_squared_inv = (1 / sin_beta_sq) * (
            h**2 / a**2 +
            (k**2 * sin_beta_sq) / b**2 +
            l**2 / c**2 -
            (2 * h * l * cos_beta) / (a * c)
        )
        
        if d_squared_inv <= 0:
            continue

        d_spacing = np.sqrt(1 / d_squared_inv)

        # Bragg's Law: 2d*sin(θ) = λ  =>  θ = arcsin(λ / 2d)
        # Argument for arcsin must be <= 1 for a real solution
        sin_theta_arg = lambda_val / (2 * d_spacing)
        if sin_theta_arg > 1:
            continue
            
        theta_rad = np.arcsin(sin_theta_arg)
        calculated_2theta = np.degrees(2 * theta_rad)

        # Format (h k l) string for nice printing
        hkl_str = f"({h} {k} {l})"
        
        # Print one row of the results table
        row = (
            f"{approx_2theta:<30.1f} "
            f"{hkl_str:<25} "
            f"{d_spacing:<30.4f} "
            f"{calculated_2theta:<20.2f}"
        )
        print(row)

# Execute the function to generate and print the table
generate_cuo_xrd_table()