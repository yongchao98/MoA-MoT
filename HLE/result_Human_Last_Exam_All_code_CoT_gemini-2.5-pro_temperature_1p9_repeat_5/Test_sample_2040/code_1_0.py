def check_spdc_potential(material_name, is_centrosymmetric, notes=""):
    """
    Analyzes the potential for a material to exhibit Spontaneous Parametric Down-Conversion (SPDC).

    Args:
        material_name (str): The name of the material.
        is_centrosymmetric (bool): True if the material has inversion symmetry, False otherwise.
        notes (str): Additional context or explanation.
    """
    print(f"Analyzing: {material_name}")
    print("-------------------------------------------------")
    print("Principle: Spontaneous Parametric Down-Conversion (SPDC) is a second-order nonlinear optical process.")
    print("This process depends on the material's second-order susceptibility, χ⁽²⁾.")
    print("A fundamental requirement for a non-zero bulk χ⁽²⁾ is that the material must lack a center of inversion (i.e., it must be non-centrosymmetric).")
    print("\nChecking material properties...")
    print(f"Material: {material_name}")
    print(f"Is centrosymmetric? {is_centrosymmetric}")

    if is_centrosymmetric:
        print("\nConclusion: The material is centrosymmetric.")
        print("Therefore, its bulk second-order susceptibility (χ⁽²⁾) is zero.")
        print(f"As a result, free-standing {material_name} is NOT expected to exhibit spontaneous parametric downconversion.")
    else:
        print("\nConclusion: The material is non-centrosymmetric.")
        print("Therefore, its bulk second-order susceptibility (χ⁽²⁾) can be non-zero.")
        print(f"As a result, {material_name} would be expected to exhibit spontaneous parametric downconversion.")

    if notes:
        print("\nAdditional Notes:")
        print(notes)
    print("-------------------------------------------------")

if __name__ == "__main__":
    # Properties of common, stable phases of borophene
    borophene_notes = ("The most stable and commonly studied phases of borophene (e.g., β₁₂ and χ₃) have crystal structures with a center of inversion. "
                       "While symmetry-breaking effects like strain or defects could potentially induce a weak second-order response, an ideal free-standing sheet would not exhibit SPDC.")

    check_spdc_potential(
        material_name="Boron Nanosheet (Borophene)",
        is_centrosymmetric=True,
        notes=borophene_notes
    )