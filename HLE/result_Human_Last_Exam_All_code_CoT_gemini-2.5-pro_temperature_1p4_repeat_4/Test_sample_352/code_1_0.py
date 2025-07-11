def suggest_mutagenesis():
    """
    Suggests a site-directed mutagenesis plan to relieve the negative charge
    from a specified protein patch.
    """

    # Define the original amino acids and their positions
    original_patch = {
        47: 'S',
        48: 'E',
        49: 'E',
        50: 'D'
    }

    # Define the proposed replacement amino acid
    replacement_aa = 'A'
    replacement_name = 'Alanine'

    # Explain the reasoning
    print("Designing a site-directed mutagenesis experiment to neutralize the negatively charged patch (positions 47-50).")
    print("The goal is to replace the original amino acids with a small, neutral amino acid to test the role of the negative charge in autoinhibition.")
    print(f"The best replacement candidate is {replacement_name} ({replacement_aa}) because it is uncharged and minimally disruptive.\n")
    
    print("The proposed change replaces the original S-E-E-D sequence with A-A-A-A.")
    print("This mutation will eliminate the negative charges from E and D, and also prevent the potential negative charge from the phosphorylation of S.\n")

    # Build the 'equation' parts
    original_parts = []
    mutated_parts = []

    # Sort positions to ensure consistent order
    for position in sorted(original_patch.keys()):
        original_aa = original_patch[position]
        original_parts.append(f"{position}{original_aa}")
        mutated_parts.append(f"{position}{replacement_aa}")

    # Format the final equation string
    original_str = " + ".join(original_parts)
    mutated_str = " + ".join(mutated_parts)

    # Print the final mutation plan as an equation
    print("Final mutation plan:")
    print(f"{original_str}  ->  {mutated_str}")

# Execute the function to print the plan
suggest_mutagenesis()