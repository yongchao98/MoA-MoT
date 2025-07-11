import math

def solve_pioneer_probe_packing():
    """
    Solves the Pioneer probe packing problem by searching for the optimal
    box dimensions and ball arrangement.
    """
    
    # Use integer units of 0.5 cm to avoid floating point issues.
    # 1 cm = 2 units.
    R_LARGE_UNITS = 4  # 2-cm radius
    R_SMALL_UNITS = 2  # 1-cm radius
    SPACING_LARGE_UNITS = 2 * R_LARGE_UNITS
    
    E_LARGE = 20  # MJ
    E_SMALL = 1   # MJ
    
    MAX_SA_CM2 = 1050
    
    best_energy = 0
    best_config = None

    # Iterate through possible dimensions in 0.5 cm steps.
    # We assume L >= W >= H to avoid redundant permutations.
    # Min dimension to fit one large ball is 2*radius = 4cm.
    # A rough upper bound from a cube: 6*L^2=1050 -> L~13. Let's search wider.
    for l_cm in [i * 0.5 for i in range(8, 41)]:  # L from 4.0cm to 20.0cm
        l = int(l_cm * 2)
        for w_cm in [i * 0.5 for i in range(8, int(l_cm*2) + 1)]: # W from 4.0cm to L
            w = int(w_cm * 2)
            
            # Estimate min H to maybe satisfy SA, to prune search
            if 2*l*w / 4.0 > MAX_SA_CM2:
                continue

            for h_cm in [i * 0.5 for i in range(8, int(w_cm*2) + 1)]: # H from 4.0cm to W
                h = int(h_cm * 2)

                sa_cm2 = 2 * (l*w + l*h + w*h) / 4.0
                if sa_cm2 > MAX_SA_CM2:
                    continue

                # --- Calculate packing for this box ---
                
                # 1. Pack large balls in a simple cubic lattice
                if l < 2*R_LARGE_UNITS or w < 2*R_LARGE_UNITS or h < 2*R_LARGE_UNITS:
                    n_large = 0
                    nx_l, ny_l, nz_l = 0, 0, 0
                else:
                    nx_l = math.floor((l - 2*R_LARGE_UNITS) / SPACING_LARGE_UNITS) + 1
                    ny_l = math.floor((w - 2*R_LARGE_UNITS) / SPACING_LARGE_UNITS) + 1
                    nz_l = math.floor((h - 2*R_LARGE_UNITS) / SPACING_LARGE_UNITS) + 1
                    n_large = nx_l * ny_l * nz_l

                # 2. Pack small balls in the main interstitial sites of the large ball lattice
                n_small = 0
                if n_large > 0 and nx_l > 1 and ny_l > 1 and nz_l > 1:
                     # Check if small balls in these voids don't collide with large ones.
                     # Dist from interstitial center to large ball center is sqrt(3)*R_LARGE = 1.732 * 2 = 3.46 cm.
                     # Required distance is r_large + r_small = 2 + 1 = 3 cm. It fits.
                     n_small = (nx_l - 1) * (ny_l - 1) * (nz_l - 1)

                current_energy = n_large * E_LARGE + n_small * E_SMALL

                if current_energy > best_energy:
                    best_energy = current_energy
                    best_config = {
                        "shape": "box",
                        "l": l / 2.0, "w": w / 2.0, "h": h / 2.0,
                        "n_large": n_large,
                        "n_small": n_small,
                        "energy": current_energy,
                        "sa": sa_cm2
                    }
                    
    return best_config

# --- Main execution and output ---
result = solve_pioneer_probe_packing()

if result:
    container_desc = f"{result['shape']} {result['l']}x{result['w']}x{result['h']}"
    num_small_balls = result['n_small']
    num_large_balls = result['n_large']
    total_energy = result['energy']

    print(f"Optimal design found:")
    print(f"Container: A {container_desc}")
    print(f"Surface Area: {result['sa']:.2f} cm^2 (Constraint: <= 1050 cm^2)")
    print(f"\nPacked Balls:")
    print(f"- Number of 1-cm balls (a): {num_small_balls}")
    print(f"- Number of 2-cm balls (b): {num_large_balls}")
    
    print("\nTotal Potential Energy Calculation:")
    print(f"({num_large_balls} * {E_LARGE}) + ({num_small_balls} * {E_SMALL}) = {total_energy} MJ")
else:
    print("No valid configuration found.")
