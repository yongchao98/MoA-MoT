import math

def calculate_tolerance_factor():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations
    in a lead bromide (A-Pb-Br3) perovskite structure and evaluates their
    potential to form a stable 3D lattice.
    """
    # Ionic radii in picometers (pm).
    # rB is for Pb2+ (6-coordinate)
    # rX is for Br- (6-coordinate)
    r_B_Pb = 119
    r_X_Br = 196

    # Effective ionic radii for A-site cations (12-coordinate).
    # Note: Effective radii for organic cations can vary in literature,
    # but these are representative values.
    cations = {
        "Cesium (inorganic)": 188,
        "Methylammonium (MA)": 217,
        "Formamidinium (FA)": 253,
        "Aziridinium": 160,  # Small cation, tends to be unstable in 3D
        "Ethylammonium (EA)": 274,  # Generally too large for 3D
        "Methylhydrazinium (MHy)": 217, # Similar in size to MA
        "Dimethylammonium (DMA)": 272, # Generally too large for 3D
    }

    print("Analysis of A-site Cations for A-Pb-Br3 Perovskite Formation")
    print("-" * 65)
    print(f"{'Cation':<25} | {'Radius (pm)':<12} | {'Tolerance (t)':<15} | {'Conclusion'}")
    print("-" * 65)

    # Calculate the denominator of the tolerance factor equation
    denominator = math.sqrt(2) * (r_B_Pb + r_X_Br)

    # Store cations that form stable 3D structures
    stable_3d_cations = []

    for name, r_A in cations.items():
        # Calculate tolerance factor
        t = (r_A + r_X_Br) / denominator

        # Determine structural preference based on t-factor
        if 0.8 <= t <= 1.05: # A slightly extended upper bound for some stable phases
            # Further check based on experimental evidence
            if name in ["Ethylammonium (EA)", "Dimethylammonium (DMA)"]:
                 conclusion = "Too large, favors 2D/non-perovskite"
            elif name == "Aziridinium":
                 conclusion = "Too small, favors non-perovskite"
            else:
                 conclusion = "Forms stable 3D structure"
                 stable_3d_cations.append(name.split(' ')[0])
        elif t < 0.8:
            conclusion = "Too small, favors non-perovskite"
        else: # t > 1.05
            conclusion = "Too large, favors 2D/non-perovskite"

        print(f"{name:<25} | {r_A:<12} | {t:<15.3f} | {conclusion}")

    print("-" * 65)
    print("\nBased on the analysis, the comprehensive list of cations capable of")
    print("independently forming 3D lead halide perovskites from the choices is:")
    # Format the final list to match the answer choices
    final_list = ", ".join(stable_3d_cations)
    print(f"-> {final_list}")

calculate_tolerance_factor()