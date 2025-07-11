def check_spdc_in_boron_nanosheets():
    """
    Explains whether free-standing boron nanosheets are expected to exhibit
    spontaneous parametric down-conversion (SPDC).
    """
    # Define the core requirement for SPDC
    spdc_requirement = "a non-centrosymmetric crystal structure (lacking a center of inversion)"
    
    # Define the property of centrosymmetric materials
    centrosymmetric_property = "the second-order nonlinear susceptibility (χ⁽²⁾) is zero"
    
    # Describe the structure of common boron nanosheets (borophene)
    borophene_structure = "centrosymmetric"
    
    # Formulate the conclusion
    conclusion = (
        f"1. Spontaneous Parametric Down-Conversion (SPDC) requires {spdc_requirement}.\n"
        f"2. In materials with a center of inversion, {centrosymmetric_property}, which forbids SPDC.\n"
        f"3. Common, stable phases of boron nanosheets have a {borophene_structure} crystal structure.\n"
        f"4. Conclusion: Therefore, free-standing boron nanosheets are generally NOT expected to exhibit SPDC."
    )
    
    print(conclusion)

# Execute the function to print the explanation.
check_spdc_in_boron_nanosheets()