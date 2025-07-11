def analyze_spdc_in_boron_nanosheets():
    """
    Analyzes whether boron nanosheets are expected to exhibit SPDC
    by checking the physical requirements.
    """
    print("Analyzing the potential for Spontaneous Parametric Down-Conversion (SPDC) in Boron Nanosheets...")
    print("="*80)

    # Step 1: Explain the primary requirement for SPDC.
    print("\nStep 1: The Requirement for Second-Order Nonlinearity")
    print("Spontaneous Parametric Down-Conversion (SPDC) is a second-order nonlinear optical process.")
    print("For a material to exhibit SPDC, it must have a significant, non-zero second-order nonlinear susceptibility, denoted as χ^(2) (chi-2).")
    print("-" * 80)

    # Step 2: Connect the optical requirement to crystal symmetry.
    print("\nStep 2: The Role of Crystal Symmetry")
    print("A fundamental principle in optics is that a material can only possess a non-zero χ^(2) if its crystal structure lacks a center of inversion.")
    print("Materials with a center of inversion are called 'centrosymmetric', and they have a χ^(2) of zero in the electric-dipole approximation.")
    print("-" * 80)

    # Step 3: Assess the symmetry of boron nanosheets.
    print("\nStep 3: Crystal Structure of Boron Nanosheets (Borophene)")
    print("Common, stable phases of ideal, free-standing boron nanosheets (e.g., β₁₂ and χ₃ phases) are known to be centrosymmetric.")
    print("This means they possess a center of inversion.")
    print("-" * 80)

    # Step 4: Draw the conclusion using the governing equation.
    print("\nStep 4: Conclusion and Governing Equation")
    print("The second-order nonlinear polarization P^(2), which drives effects like SPDC, is given by:")
    print("\n\tP^(2) = ε₀ * χ^(2) * E²\n")
    print("Where:")
    print("\t- P^(2) is the second-order polarization.")
    print("\t- ε₀ is the permittivity of free space.")
    print("\t- E is the electric field of the incident light.")
    print("\t- χ^(2) is the second-order susceptibility.")
    
    # Define the value for the material in question.
    chi_2_value_for_borophene = 0
    
    print("\nFor ideal boron nanosheets, since the material is centrosymmetric:")
    print(f"\tThe value of χ^(2) is {chi_2_value_for_borophene}.")
    
    print("\nPlugging this into the equation:")
    final_polarization = 0
    print(f"\tP^(2) = ε₀ * {chi_2_value_for_borophene} * E² = {final_polarization}")
    
    print("\nBecause the second-order polarization P^(2) is zero, ideal, free-standing boron nanosheets are not expected to exhibit SPDC.")
    print("="*80)

    print("\nNOTE: This conclusion can change if the material's symmetry is broken by external factors like applied strain, structural defects, or strong electric fields.")


# Run the analysis
analyze_spdc_in_boron_nanosheets()