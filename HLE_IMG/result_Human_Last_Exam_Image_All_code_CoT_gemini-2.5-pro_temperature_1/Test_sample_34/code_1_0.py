import pandas as pd

def solve_kinship_diagram():
    """
    Interprets a Lévi-Strauss kinship diagram by comparing its relational
    pattern to ethnographic data.
    """
    # Step 1 & 2: Define the kinship data from Lévi-Strauss's analysis.
    # +1 represents a familiar (+) relationship.
    # -1 represents a formal/distant (-) relationship.
    # The relationships are Brother/Sister (B/S), Husband/Wife (H/W),
    # Father/Son (F/S), and Mother's Brother/Sister's Son (MB/ZS).
    data = {
        'Society': ['Trobriand', 'Cherkess', 'Tonga', 'Siuoi'],
        'Lineage': ['Matrilineal', 'Patrilineal', 'Patrilineal', 'Matrilineal'],
        'B/S': [-1, 1, -1, 1],
        'H/W': [1, -1, 1, -1],
        'F/S': [1, -1, -1, 1],
        'MB/ZS': [-1, 1, 1, -1]
    }
    kinship_df = pd.DataFrame(data).set_index('Society')

    # Step 3: Define the pattern from the diagram.
    # Father/Son (F/S) line has a '-' sign.
    # Mother's Brother/Sister's Son (MB/ZS) line has a '+' sign.
    diagram_fs = -1
    diagram_mbzs = 1

    print("Interpreting the Lévi-Strauss Kinship Diagram:")
    print("The diagram shows four key relationships:")
    print("  - Brother/Sister: - (Formal)")
    print("  - Husband/Wife: + (Familiar)")
    print("  - Father/Son: - (Formal)")
    print("  - Mother's Brother/Sister's Son: + (Familiar)\n")
    
    print("The core of the analysis focuses on the 'descent axis' relationships.")
    print(f"The Father/Son relationship in the diagram is formal, which we represent as the number: {diagram_fs}")
    print(f"The Mother's Brother/Sister's Son relationship is familiar, represented as the number: {diagram_mbzs}\n")

    print("Searching for societies that match this descent axis pattern...")
    
    # Step 4: Filter the societies based on the diagram's descent axis.
    matching_societies = kinship_df[
        (kinship_df['F/S'] == diagram_fs) &
        (kinship_df['MB/ZS'] == diagram_mbzs)
    ]

    if not matching_societies.empty:
        for society, row in matching_societies.iterrows():
            print(f"- Match Found: {society} ({row['Lineage']})")
            print(f"  - Father/Son relationship: {row['F/S']}")
            print(f"  - Mother's Brother/Sister's Son relationship: {row['MB/ZS']}")
        
        print("\nConclusion: The societies of Tonga and Cherkess both exhibit the descent structure shown in the diagram.")
        print("This corresponds to answer choice D.")
    else:
        print("No matching societies found based on the descent axis.")

solve_kinship_diagram()