import pandas as pd

def solve_kinship_diagram():
    """
    Analyzes a Lévi-Strauss kinship diagram to identify the societies it represents.
    """
    print("Step 1: Deconstructing the Lévi-Strauss kinship diagram.")
    print("The diagram uses the following symbols:")
    print("  Δ : Male")
    print("  o : Female")
    print("  = : Marriage bond (Husband-Wife)")
    print("  Horizontal line: Sibling bond (Brother-Sister)")
    print("  Diagonal line (downwards): Descent bond (Parent-Son)")
    print("  + : Familiar / Affectionate relationship")
    print("  - : Formal / Authoritarian / Antagonistic relationship\n")

    print("Step 2: Interpreting the relationships and attitudes from the diagram.")
    
    # Diagram analysis
    # Let's label the individuals: Left Δ (Father), Center o (Wife/Mother), Right Δ (Mother's Brother), Bottom Δ (Son)
    # Husband-Wife: The '=' line has a '-' sign.
    husband_wife_attitude = "-"
    # Father-Son: The line from the Father (left Δ) to the Son (bottom Δ) has a '+' sign.
    father_son_attitude = "+"
    # Brother-Sister: The line connecting the sister (o) and brother (right Δ) has '+' and '-' signs.
    # In Lévi-Strauss's theory, the B-S relationship is opposite to the H-W relationship.
    # Since H-W is '-', B-S must be '+'.
    brother_sister_attitude = "+"
    # Mother's Brother - Sister's Son: This relationship's attitude is opposite to the Father-Son relationship.
    # Since F-S is '+', MB-SS must be '-'.
    mb_ss_attitude = "-"

    print(f"From the diagram, we deduce the following relationship attitudes:")
    print(f"  - Husband-Wife: {husband_wife_attitude} (Formal/Distant)")
    print(f"  - Brother-Sister: {brother_sister_attitude} (Familiar/Close)")
    print(f"  - Father-Son: {father_son_attitude} (Familiar/Close)")
    print(f"  - Mother's Brother - Sister's Son: {mb_ss_attitude} (Formal/Authoritarian)\n")

    print("Step 3: Determining the descent system represented.")
    print("The key indicator is the axis of authority between the son and the males of the previous generation.")
    print(f"Here, the Father-Son relationship is familiar ({father_son_attitude}), while the Mother's Brother - Sister's Son relationship is one of authority ({mb_ss_attitude}).")
    print("Authority, discipline, and inheritance rights pass through the mother's line (from her brother).")
    print("This is the classic structure of a MATRILINEAL system.\n")

    print("Step 4: Evaluating the answer choices based on known ethnographic data.")
    societies = {
        'Trobriand': {'type': 'Matrilineal', 'Fa-So': '+', 'MB-SS': '-'},
        'Siuoi': {'type': 'Matrilineal', 'Fa-So': '+', 'MB-SS': '-'},
        'Lake Kubutu': {'type': 'Patrilineal', 'Fa-So': '-', 'MB-SS': '+'},
        'Tonga': {'type': 'Patrilineal', 'Fa-So': '-', 'MB-SS': '+'},
        'Cherkess': {'type': 'Patrilineal', 'Fa-So': '-', 'MB-SS': '+'}
    }
    
    # Using pandas for a clean table display
    df = pd.DataFrame(societies).T
    df.index.name = 'Society'
    df.columns = ['System Type', 'Father-Son Attitude', "Mother's Brother-Son Attitude"]
    print(df)
    print("\nThe diagram represents a system with a (+) Father-Son and a (-) Mother's Brother-Son relationship.")
    print("Based on the table, only Matrilineal systems like the Trobriand and Siuoi fit this pattern.\n")

    print("Step 5: Selecting the correct answer choice.")
    print("The choices are:")
    print("A. Trobriand-matrilineal and Siuoi-matrilineal")
    print("B. Siuoi-matrilineal and Lake Kubutu-patrilineal")
    print("C. Lake Kubutu-patrilineal and Tonga-patrilineal")
    print("D. Tonga-patrilineal and Cherkess-patrilineal")
    print("E. Cherkess-patrilineal and Trobiand-matrilineal")
    
    print("\nOnly choice A contains two societies whose matrilineal structure matches the diagram.")

if __name__ == '__main__':
    solve_kinship_diagram()