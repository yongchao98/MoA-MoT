import sys

def check_spdc_capability(material_name, is_centrosymmetric):
    """
    Checks if a material is expected to exhibit Spontaneous Parametric
    Down-Conversion (SPDC) based on its crystal symmetry.

    The governing principle is that SPDC is a second-order nonlinear
    optical effect. Its strength is dependent on the second-order
    susceptibility, chi^(2). In materials with a center of inversion
    (centrosymmetric materials), chi^(2) is required by symmetry to be zero.
    """
    # Use UTF-8 encoding to ensure special characters print correctly
    sys.stdout.reconfigure(encoding='utf-8')

    # Define unicode characters for the formula to be printed
    chi = '\u03C7'
    sup_l_paren = '\u207D'
    sup_2 = '\u00B2'
    sup_r_paren = '\u207E'
    chi_2_symbol = f"{chi}{sup_l_paren}{sup_2}{sup_r_paren}"

    print(f"Analysis for: {material_name}")
    print(f"Material Property (Is Centrosymmetric): {is_centrosymmetric}")
    print("-" * 40)

    # In a simplified model, the effect's strength is proportional to chi^(2).
    # Equation: Effect_Strength = C * chi^(2), where C is a constant.
    # The key number in this equation is the value of chi^(2).
    if is_centrosymmetric:
        chi_2_value = 0
    else:
        # For a non-centrosymmetric material, chi^(2) would be non-zero.
        # We use a non-zero placeholder for demonstration purposes.
        chi_2_value = 1.0  # Using an arbitrary value in placeholder units

    print(f"The simplified governing equation is: SPDC_Effect \u221D {chi_2_symbol}")
    print(f"For this material, the key number in the equation is:")
    print(f"Value of {chi_2_symbol} = {chi_2_value}")

    if chi_2_value == 0:
        print(f"\nConclusion: Because {chi_2_symbol} is 0, the material is not expected to exhibit SPDC.")
    else:
        print(f"\nConclusion: Because {chi_2_symbol} is not 0, the material could exhibit SPDC.")


# Most stable, free-standing phases of boron nanosheets are centrosymmetric.
material_to_check = {
    "name": "Free-standing Boron Nanosheet (common phases)",
    "is_centrosymmetric": True
}

# Run the analysis for the boron nanosheet.
check_spdc_capability(material_to_check["name"], material_to_check["is_centrosymmetric"])