def check_spdc_capability(material_name, is_centrosymmetric):
    """
    Checks if a material is expected to exhibit Spontaneous Parametric Down-Conversion (SPDC)
    based on its crystal symmetry.
    """

    print("Analyzing the potential for Spontaneous Parametric Down-Conversion (SPDC) in:", material_name)
    print("-" * 70)

    # Step 1: Explain the fundamental requirement for SPDC.
    print("Step 1: SPDC is a second-order nonlinear optical effect.")
    print("         Its occurrence depends on a material property called the second-order")
    print("         nonlinear susceptibility, denoted as chi_2 (χ⁽²⁾).")

    # Step 2: Explain the role of material symmetry.
    print("\nStep 2: A material's crystal symmetry dictates whether chi_2 is zero or non-zero.")
    print("         For materials that have a center of inversion (are centrosymmetric),")
    print("         chi_2 is required by symmetry to be zero.")
    print("         Therefore, only non-centrosymmetric materials can exhibit SPDC.")

    # Step 3: Apply this principle to the specific material.
    print(f"\nStep 3: Evaluating the crystal structure of {material_name}.")
    if is_centrosymmetric:
        print(f"         Research shows that the common, stable phases of free-standing")
        print(f"         {material_name} (borophene) are centrosymmetric.")
        chi_2 = 0
    else:
        # This case is not applicable for borophene but is included for completeness.
        print(f"         {material_name} is found to be non-centrosymmetric.")
        chi_2 = "non-zero" 

    # Step 4: Formulate the conclusion based on the derived value of chi_2.
    print("\nStep 4: Based on the symmetry, the value for the key parameter is determined.")
    print(f"         The equation for the bulk second-order susceptibility is: chi_2 = {chi_2}")

    if chi_2 == 0:
        conclusion = (
            f"Conclusion: No. Bulk, free-standing {material_name} are NOT expected to exhibit\n"
            "            spontaneous parametric downconversion. This is because they are\n"
            "            centrosymmetric, which forces their second-order susceptibility (chi_2)\n"
            "            to be zero, forbidding the effect."
        )
    else:
        conclusion = (
            f"Conclusion: Yes. {material_name} would be expected to exhibit SPDC because it is\n"
            "            non-centrosymmetric, allowing for a non-zero chi_2."
        )

    print("\n" + "="*70)
    print(conclusion)
    print("="*70)

# Define the properties for free-standing boron nanosheets
material = "boron nanosheets"
# Set its known symmetry property.
symmetry_is_centrosymmetric = True

# Run the analysis.
check_spdc_capability(material, symmetry_is_centrosymmetric)