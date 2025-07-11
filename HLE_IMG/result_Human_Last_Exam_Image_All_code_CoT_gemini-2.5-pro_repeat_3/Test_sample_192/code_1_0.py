import numpy as np

def index_cuo_xrd():
    """
    Calculates and prints a table of indexed XRD peaks for monoclinic CuO.
    """
    print("Indexing XRD Peaks for CuO Nanoparticles")
    print("=====================================================")
    print("This script identifies the Miller indices (h k l) for the prominent peaks in")
    print("the provided XRD spectrum of CuO. The calculations are based on the monoclinic")
    print("crystal structure of CuO (Tenorite) and Bragg's Law.")
    print("\n--- Parameters Used ---")

    # Lattice parameters for CuO (Tenorite, monoclinic, JCPDS 05-0661)
    a = 4.683  # in Angstrom
    b = 3.425  # in Angstrom
    c = 5.129  # in Angstrom
    beta_deg = 99.47  # in degrees
    
    # X-ray wavelength for Cu K-alpha radiation
    lambda_ = 1.5406  # in Angstrom

    print(f"Lattice Parameters: a = {a} Å, b = {b} Å, c = {c} Å, β = {beta_deg}°")
    print(f"X-ray Wavelength (λ): {lambda_} Å")
    print("-----------------------------------------------------\n")

    def calculate_2theta(h, k, l):
        """Calculates the 2-theta angle for a given (h,k,l) plane."""
        beta_rad = np.deg2rad(beta_deg)
        
        # Equation for interplanar spacing (d) in a monoclinic system
        # 1/d^2 = (1/sin^2(β)) * [h^2/a^2 + k^2*sin^2(β)/b^2 + l^2/c^2 - 2*h*l*cos(β)/(a*c)]
        sin_beta_sq = np.sin(beta_rad)**2
        cos_beta = np.cos(beta_rad)
        
        d_inv_sq = (1 / sin_beta_sq) * (
            (h**2 / a**2) + 
            (k**2 * sin_beta_sq / b**2) + 
            (l**2 / c**2) - 
            (2 * h * l * cos_beta / (a * c))
        )
        
        if d_inv_sq <= 0:
            return None
        
        d_spacing = np.sqrt(1 / d_inv_sq)
        
        # Bragg's Law: 2d*sin(theta) = λ
        # sin(theta) = λ / (2*d)
        sin_theta = lambda_ / (2 * d_spacing)
        
        if sin_theta > 1:
            return None
            
        theta_rad = np.arcsin(sin_theta)
        two_theta_deg = np.rad2deg(2 * theta_rad)
        return two_theta_deg

    # Prominent peaks identified from the plot and their known (h,k,l) indices from literature
    peak_data = [
        (32.5, (1, 1, 0)),
        (35.5, (-1, 1, 1)),
        (38.7, (1, 1, 1)),
        (48.7, (2, 0, -2)),
        (53.5, (0, 2, 0)),
        (58.3, (2, 0, 2)),
    ]

    # --- Generate and Print the Final Table ---
    print("--- Peak Indexing Table ---")
    header = f"{'Experimental 2θ (°)':<25} {'(h k l) Index':<18} {'Calculated 2θ (°)':<20}"
    print(header)
    print("-" * len(header))

    for exp_2theta, hkl in peak_data:
        h, k, l = hkl
        calc_2theta = calculate_2theta(h, k, l)
        
        # Format the output line for the table
        hkl_str = f"({h} {k} {l})"
        line = f"{exp_2theta:<25.1f} {hkl_str:<18} {calc_2theta:<20.2f}"
        print(line)

    print("-" * len(header))

# Execute the function
index_cuo_xrd()