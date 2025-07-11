def compare_biopolymers():
    """
    Prints a detailed explanation of why polynucleotides are not,
    structurally, polysaccharides.
    """
    direct_answer = "No, polynucleotides are not structurally polysaccharides."
    
    # Header
    print("The Question: Are homopolynucleotides, structurally, polysaccharides?\n")
    print(f"The direct answer is: {direct_answer}\n")
    print("While both are large biological polymers, their chemical structures are fundamentally different.")
    print("Here is a breakdown of the key structural differences:\n")

    # Data for comparison table
    comparison_data = [
        ("Feature", "Polynucleotide (e.g., DNA, RNA)", "Polysaccharide (e.g., Starch, Cellulose)"),
        ("="*15, "="*35, "="*45),
        ("Monomer", "Nucleotide", "Monosaccharide (simple sugar)"),
        ("Monomer Parts", "1. Phosphate Group\n\t\t 2. 5-Carbon Sugar\n\t\t 3. Nitrogenous Base", "A single sugar unit (e.g., glucose)"),
        ("Backbone", "Sugar-Phosphate Backbone", "Chain of sugar units"),
        ("Linkage Bond", "Phosphodiester Bond", "Glycosidic Bond"),
        ("Atoms Present", "C, H, O, N, P", "C, H, O")
    ]
    
    # Print the formatted table
    print(f"{comparison_data[0][0]:<15} | {comparison_data[0][1]:<35} | {comparison_data[0][2]:<45}")
    print(f"{comparison_data[1][0]:<15} | {comparison_data[1][1]:<35} | {comparison_data[1][2]:<45}")

    # Print data rows, handling multiline entries
    for feature, poly_n_desc, poly_s_desc in comparison_data[2:]:
        # Split multiline descriptions for aligned printing
        poly_n_lines = poly_n_desc.split('\n')
        poly_s_lines = poly_s_desc.split('\n')
        
        max_lines = max(len(poly_n_lines), len(poly_s_lines))
        
        for i in range(max_lines):
            feature_part = feature if i == 0 else ""
            poly_n_part = poly_n_lines[i] if i < len(poly_n_lines) else ""
            poly_s_part = poly_s_lines[i] if i < len(poly_s_lines) else ""
            print(f"{feature_part:<15} | {poly_n_part:<35} | {poly_s_part:<45}")

    # Conclusion
    print("\n" + "="*80)
    print("\nConclusion:")
    print("The only structural similarity is that both contain sugar molecules. However, in polynucleotides,")
    print("the sugar is just one part of a more complex nucleotide monomer. The presence of phosphate groups")
    print("in the backbone and nitrogenous bases as side groups makes polynucleotides a distinct chemical")
    print("class from polysaccharides, which are composed solely of sugar units.")

# Run the explanation
compare_biopolymers()
