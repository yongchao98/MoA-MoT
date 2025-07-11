def explain_spdc_in_borophene():
    """
    Explains whether free-standing boron nanosheets are expected to exhibit
    spontaneous parametric downconversion (SPDC).
    """

    print("--- Step 1: Requirements for Spontaneous Parametric Downconversion (SPDC) ---")
    print("SPDC is a second-order nonlinear optical process where a photon splits into two photons of lower energy.")
    print("A fundamental requirement for a material to facilitate this process is a non-zero second-order nonlinear susceptibility, denoted as χ⁽²⁾.")
    print("For χ⁽²⁾ to be non-zero, the material's crystal structure must be non-centrosymmetric (i.e., lack a center of inversion).")
    print("\n")

    print("--- Step 2: Crystal Structure of Boron Nanosheets (Borophene) ---")
    print("Free-standing boron nanosheets are known to have multiple stable structures (polymorphism).")
    print("The most commonly synthesized and studied phases, such as the β12 and χ3 phases, have a pmmn space group.")
    print("A material with the pmmn space group has a crystal structure that is centrosymmetric (it possesses a center of inversion).")
    print("\n")

    print("--- Step 3: Conclusion ---")
    print("Since the common, stable phases of free-standing boron nanosheets are centrosymmetric, their intrinsic χ⁽²⁾ is zero.")
    print("Therefore, ideal, free-standing boron nanosheets are NOT expected to exhibit spontaneous parametric downconversion.")
    print("\nCaveat: While ideal nanosheets are not expected to show this effect, the introduction of strain, defects, or interactions with a substrate could potentially break the inversion symmetry, making a weak second-order response possible in non-ideal conditions.")

if __name__ == '__main__':
    explain_spdc_in_borophene()