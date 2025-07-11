import numpy as np

def solve_chemistry_puzzle():
    """
    Solves the user's puzzle by performing chemical graph theory calculations on Phosphine (PH3).
    
    The choice of PH3 is based on the following reasoning:
    1. The puzzle's clues point towards the element Phosphorus.
    2. The 'Mass-Weighted Moran's I' calculation requires a hetero-atomic molecule
       to be well-defined and non-trivial.
    3. Phosphine (PH3) is the simplest, common compound of Phosphorus, making it the most
       logical candidate for a puzzle of this nature.
    """

    # --- Step 1: Define molecule properties and geometry for Phosphine (PH3) ---
    
    # Atomic properties from IUPAC
    Z_P, m_P = 15, 30.9737620
    Z_H, m_H = 1, 1.0080

    # Molecule structure data
    d_ph = 1.419  # P-H bond length in Angstroms
    angle_hph_deg = 93.5  # H-P-H bond angle in degrees
    
    # Create a list of atoms with their properties
    # Atom 0: P, Atoms 1-3: H
    atoms_Z = [Z_P, Z_H, Z_H, Z_H]
    atoms_m = [m_P, m_H, m_H, m_H]
    
    # Calculate 3D coordinates for a trigonal pyramid
    # Based on P-H bond length and H-P-H angle
    angle_hph_rad = np.deg2rad(angle_hph_deg)
    
    # Using formulas to derive pyramid height (h) and base triangle radius (r)
    # cos(angle) = (-r^2/2 + h^2) / (r^2+h^2)  and d^2 = r^2 + h^2
    # Solving these gives h and r
    cos_angle = np.cos(angle_hph_rad)
    # Ratio h^2/r^2 from the angle formula
    h2_r2_ratio = (cos_angle + 0.5) / (1 - cos_angle)
    # Find r^2 from d^2 = r^2 + h^2 = r^2(1 + h^2/r^2)
    r2 = d_ph**2 / (1 + h2_r2_ratio)
    r = np.sqrt(r2)
    h = np.sqrt(d_ph**2 - r2)
    
    # Set up coordinates
    # P is on the z-axis, H's are in the xy-plane forming an equilateral triangle
    coords = np.array([
        [0.0, 0.0, h],                                 # P atom 0
        [r, 0.0, 0.0],                                 # H atom 1
        [r * np.cos(2*np.pi/3), r * np.sin(2*np.pi/3), 0.0],  # H atom 2
        [r * np.cos(4*np.pi/3), r * np.sin(4*np.pi/3), 0.0]   # H atom 3
    ])

    # --- Step 2: Calculate Mass-Weighted Barysz Graph Energy ---

    # Adjacency matrix (connectivity C_ij): P is bonded to 3 H's
    C = np.array([
        [0, 1, 1, 1],
        [1, 0, 0, 0],
        [1, 0, 0, 0],
        [1, 0, 0, 0]
    ])

    # Distance matrix d_ij
    d = np.zeros((4, 4))
    for i in range(4):
        for j in range(i + 1, 4):
            dist = np.linalg.norm(coords[i] - coords[j])
            d[i, j] = d[j, i] = dist

    # Barysz Matrix B_ij = C_ij * sqrt(Z_i*Z_j / (d_ij * m_i))
    B = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            if C[i, j] == 1:
                numerator = atoms_Z[i] * atoms_Z[j]
                denominator = d[i, j] * atoms_m[i]
                B[i, j] = np.sqrt(numerator / denominator)

    # Eigenvalues and Energy (sum of absolute values)
    eigenvalues = np.linalg.eigvals(B)
    barysz_energy = np.sum(np.abs(eigenvalues))

    # --- Step 3: Calculate Min/Max Mass-Weighted Moran's I ---
    
    # The variable for Moran's I is atomic mass
    x = np.array(atoms_m)
    x_bar = np.mean(x)
    n = len(x)

    # Weights w_ij = C_ij (1 if bonded, 0 otherwise)
    w = C
    
    # Calculate local Moran's I for each atom
    # I_i = (z_i / var_z) * sum_j(w_ij * z_j) where z_i = x_i - x_bar
    z = x - x_bar
    # Global variance of x
    var_x = np.sum(z**2) / n
    
    local_moran_I = []
    for i in range(n):
        # The sum over neighbors
        spatial_lag = np.sum(w[i] * z)
        # Standard deviation (or variance) term can be written in multiple ways, 
        # a common local form is I_i = z_i * spatial_lag
        # For cross-comparison, we use a normalized form z_i / m2 where m2=var
        m2 = sum(z**2)
        local_I = (n*z[i] / m2) * spatial_lag
        local_moran_I.append(local_I)
    
    # Using a common formal definition of local Moran's I for comparison
    # I_i = z_i * sum(w_ij * z_j)
    # This leads to min/max that is not scaled as per the global statistic.
    # To keep this self-contained, let's recalculate based on the simpler z_i*lag
    z_scores = (x-x.mean())/x.std()
    moran_local_vals = []
    for i in range(n):
        spatial_lag = np.sum([w[i,j] * z_scores[j] for j in range(n)])
        moran_local_vals.append(z_scores[i] * spatial_lag)

    moran_min = min(moran_local_vals)
    moran_max = max(moran_local_vals)

    # --- Step 4: Compute the Final Product ---
    final_product = barysz_energy * moran_min * moran_max
    
    # Print the equation with final values
    print(f"{barysz_energy} * {moran_min} * {moran_max} = {final_product}")

    # Final answer in the requested format
    print(f"<<<{final_product:.1f}>>>")

solve_chemistry_puzzle()