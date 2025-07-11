def check_spdc_potential(material_name, crystal_symmetry):
    """
    Analyzes if a material is a candidate for Spontaneous Parametric Down-Conversion (SPDC)
    based on its crystal symmetry.
    """
    print(f"Analyzing the potential for {material_name} to exhibit SPDC...\n")

    # Step 1: State the physical requirement for SPDC.
    print("Step 1: Requirement for Spontaneous Parametric Down-Conversion (SPDC)")
    print("="*60)
    print("SPDC is a second-order nonlinear optical process. A material can only exhibit")
    print("this effect if its crystal structure is non-centrosymmetric (lacks a center of inversion).")
    print("If a material is centrosymmetric, its second-order susceptibility (χ⁽²⁾) is zero.\n")

    # Step 2: State the relevant property of the material.
    print(f"Step 2: Crystal Symmetry of {material_name}")
    print("="*60)
    print(f"The most stable and commonly synthesized forms of free-standing {material_name}")
    print(f"(e.g., β₁₂ and χ₃ borophene) are known to be {crystal_symmetry}.\n")

    # Step 3: Conclude based on the logic.
    print("Step 3: Conclusion")
    print("="*60)
    if crystal_symmetry == "centrosymmetric":
        print(f"Condition: Material is {crystal_symmetry}.")
        print("Result: The second-order nonlinear susceptibility (χ⁽²⁾) is zero.")
        print(f"Conclusion: Therefore, ideal free-standing {material_name} would NOT be expected to")
        print("exhibit spontaneous parametric down-conversion.")
    else:
        print(f"Condition: Material is {crystal_symmetry}.")
        print("Result: The second-order nonlinear susceptibility (χ⁽²⁾) can be non-zero.")
        print(f"Conclusion: Therefore, {material_name} would be a potential candidate for")
        print("exhibiting spontaneous parametric down-conversion.")

# --- Parameters for the specific question ---
material = "Boron Nanosheets (Borophene)"
symmetry = "centrosymmetric"

# --- Run the analysis ---
check_spdc_potential(material, symmetry)