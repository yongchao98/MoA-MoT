import collections

def analyze_kinship_diagram():
    """
    Analyzes the Lévi-Strauss kinship diagram and matches it to known systems.
    """
    
    # Step 1 & 2: Deconstruct the diagram and define the core relationships.
    # The diagram consists of a husband (Δ), wife (o), her brother (Δ), and the couple's son (Δ).
    # The signs (+ or -) indicate the nature of the relationship.
    # From visual inspection of the diagram:
    # Husband <=> Wife bond ('=') has a '-' sign.
    # Brother <=> Sister bond ('-') has a '+' sign.
    # The descent line from the marriage to the son has a '+' on the father's side.
    # The relationship between the mother's brother and the son is depicted as '-', the opposite of the father/son relationship.
    
    diagram_relationships = {
        "Husband/Wife": "-",
        "Brother/Sister": "+",
        "Father/Son": "+",
        "Mother's Brother/Sister's Son": "-"
    }

    print("--- Step 1 & 2: Decoding the Lévi-Strauss Diagram ---")
    print("The diagram represents a kinship structure with the following relationship dynamics:")
    for relation, attitude in diagram_relationships.items():
        attitude_desc = "Familiar/Informal" if attitude == "+" else "Antagonistic/Formal"
        print(f"- {relation}: {attitude_desc} ({attitude})")
    print("\n" + "="*50 + "\n")

    # Step 3: Characterize the kinship structure based on the decoded relationships.
    print("--- Step 3: Characterizing the Kinship Structure ---")
    fs_relation = diagram_relationships["Father/Son"]
    mb_zs_relation = diagram_relationships["Mother's Brother/Sister's Son"]

    print(f"The key dynamic is the relationship between paternal and maternal lines:")
    print(f" - The Father/Son relationship is familiar ({fs_relation}).")
    print(f" - The Mother's Brother/Sister's Son relationship is formal/antagonistic ({mb_zs_relation}).")
    print("\nThis pattern, where the biological father is an affectionate figure and the maternal uncle holds authority and discipline, is a classic characteristic of MATRILINEAL systems.")
    print("In such systems, descent, property, and power are inherited through the mother's line, making her brother the key male authority figure for her children.")
    print("\nConversely, in PATRILINEAL systems, we often find the inverse: a formal relationship with the father (authority) and a familiar one with the mother's brother (indulgent uncle).")
    print("\n" + "="*50 + "\n")

    # Step 4: Evaluate the societal examples from the options.
    # This data is based on anthropological ethnographies.
    known_systems = {
        "Trobriand": {"type": "matrilineal", "F/S": "+", "MB/ZS": "-"},
        "Siuoi": {"type": "matrilineal", "F/S": "+", "MB/ZS": "-"},
        "Lake Kubutu": {"type": "patrilineal", "F/S": "-", "MB/ZS": "+"},
        "Tonga": {"type": "patrilineal", "F/S": "-", "MB/ZS": "+"},
        "Cherkess": {"type": "patrilineal", "F/S": "-", "MB/ZS": "+"},
    }
    
    print("--- Step 4: Evaluating the Societal Examples ---")
    print("We will now compare the diagram's pattern with known anthropological cases:")
    
    matching_systems = []
    for society, data in known_systems.items():
        is_match = (data["F/S"] == fs_relation) and (data["MB/ZS"] == mb_zs_relation)
        print(f"\nAnalyzing: {society} ({data['type']})")
        print(f" - Father/Son relationship: {data['F/S']}")
        print(f" - Mother's Brother/Son relationship: {data['MB/ZS']}")
        print(f" - Match with diagram: {'YES' if is_match else 'NO'}")
        if is_match:
            matching_systems.append(society)

    print("\n" + "="*50 + "\n")

    # Step 5: Identify the correct match from the answer choices.
    print("--- Step 5: Final Conclusion ---")
    print("The societies that are correctly represented by the diagram are those with a matrilineal structure that matches the (+, -) pattern for Father/Son and Mother's Brother/Son relationships.")
    print(f"Our analysis shows that {', '.join(matching_systems)} fit this description.")
    print("\nLooking at the answer choices:")
    print("A. Trobriand-matrilineal and Siuoi-matrilineal")
    print("B. Siuoi-matrilineal and Lake Kubutu-patrilineal")
    print("C. Lake Kubutu-patrilineal and Tonga-patrilineal")
    print("D. Tonga-patrilineal and Cherkess-patrilineal")
    print("E. Cherkess-patrilineal and Trobiand-matrilineal")
    
    print("\nChoice A correctly identifies two systems, Trobriand and Siuoi, that match the diagram's structure.")


if __name__ == '__main__':
    analyze_kinship_diagram()