import re

def analyze_reaction():
    """
    Analyzes the specified reaction and presents the starting material,
    product, and molecular details.
    """
    # Define names and formula for the reaction components
    # The reaction is a rearrangement, so the molecular formula remains the same.
    starting_material_name = "(1S,2R,4S)-2-((S)-4-((tert-butyldimethylsilyl)oxy)cyclopent-1-en-1-yl)-7,7-dimethoxybicyclo[2.2.1]hept-5-en-2-ol"
    product_description = "A fused bicyclic ketone (structural isomer of the starting material)"
    molecular_formula = "C20H34O4Si"

    # Define the stoichiometry of the reaction equation
    reactant_coefficient = 1
    product_coefficient = 1

    # Print a summary of the transformation
    print("Reaction Analysis:")
    print("-" * 20)
    print(f"Starting Material: {starting_material_name}")
    print("Conditions: 1. KH, THF; 2. H2O/MeOH workup")
    print("Reaction Type: Anionic Oxy-Cope Rearrangement")
    print(f"Product: {product_description}")
    print("-" * 20)

    # Construct and print the final reaction equation
    final_equation = f"{reactant_coefficient} * SM ({molecular_formula}) --> {product_coefficient} * Product ({molecular_formula})"
    print("Final Equation:")
    print(f"Where SM = Starting Material and Product is its isomer.")
    print(final_equation)
    print()

    # Explicitly output each number from the equation as requested
    print("Numbers in the final equation:")
    
    # 1. Print coefficients
    print(f"Stoichiometric coefficients: {reactant_coefficient}, {product_coefficient}")

    # 2. Parse the formula and print atom counts
    atom_counts = re.findall(r'([A-Z][a-z]*)(\d*)', molecular_formula)
    elements = []
    counts = []
    for element, count in atom_counts:
        elements.append(element)
        # If no number follows an element, its count is 1.
        counts.append(int(count) if count else 1)
    
    element_str = ", ".join(elements)
    count_str = ", ".join(map(str, counts))
    print(f"Atom counts in the formula ({element_str}): {count_str}")


# Run the analysis
analyze_reaction()