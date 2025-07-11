import pandas as pd

def find_optimal_catalyst():
    """
    Simulates the process of selecting an optimal bifunctional catalyst
    for both olefin polymerization and polyolefin hydrogenolysis.
    """
    # Step 1: Model the catalyst candidates and their hypothetical performance scores.
    # Scores are on a scale of 1-10 based on known chemical principles.
    # - Polymerization: Favors stable metallocene-like structures.
    # - Hydrogenolysis: Favors highly reactive sites, often generated on acidic supports.
    data = {
        'metal': [
            'Zr', 'Zr', 'Zr',
            'Hf', 'Hf', 'Hf',
            'Ti', 'Ti', 'Ti'
        ],
        'ligand': [
            'Metallocene (Cp*2)', 'Constrained Geometry', 'Simple Alkyl (e.g., neopentyl)',
            'Metallocene (Cp*2)', 'Constrained Geometry', 'Simple Alkyl (e.g., neopentyl)',
            'Metallocene (Cp*2)', 'Constrained Geometry', 'Simple Alkyl (e.g., neopentyl)'
        ],
        'support': [
            'None', 'Silica (SiO2)', 'Sulfated Alumina (S-Al2O3)',
            'None', 'Silica (SiO2)', 'Sulfated Alumina (S-Al2O3)',
            'None', 'Silica (SiO2)', 'Sulfated Alumina (S-Al2O3)'
        ],
        # Polymerization scores
        'poly_score': [
            9, 8, 4,
            8, 9, 4,
            7, 6, 3
        ],
        # Hydrogenolysis (C-C cleavage) scores
        'hydro_score': [
            1, 3, 9,
            1, 3, 8,
            2, 4, 7
        ]
    }
    catalysts = pd.DataFrame(data)

    # Step 2: Calculate the combined score for bifunctionality.
    # A multiplicative score rewards catalysts that are good at BOTH tasks.
    catalysts['combined_score'] = catalysts['poly_score'] * catalysts['hydro_score']

    # Step 3: Identify the catalyst with the highest combined score.
    optimal_catalyst_index = catalysts['combined_score'].idxmax()
    optimal_catalyst = catalysts.loc[optimal_catalyst_index]

    # Step 4: Output the results clearly.
    print("--- Catalyst Optimization Search ---")
    print("\nSearching for a single-site catalyst effective for both polymerization and hydrogenolysis...")
    print("\n--- Optimal Catalyst Found ---")
    print(f"Metal: {optimal_catalyst['metal']}")
    print(f"Ligand: {optimal_catalyst['ligand']}")
    print(f"Support: {optimal_catalyst['support']}")
    print("-" * 30)
    print("Performance Analysis:")
    print(f"Polymerization Score: {optimal_catalyst['poly_score']}")
    print(f"Hydrogenolysis Score: {optimal_catalyst['hydro_score']}")
    print("\n--- Final Score Calculation ---")
    # Fulfilling the request to output the numbers in the final equation
    final_score = optimal_catalyst['poly_score'] * optimal_catalyst['hydro_score']
    print(f"Equation: {optimal_catalyst['poly_score']} (Polymerization) * {optimal_catalyst['hydro_score']} (Hydrogenolysis) = {final_score}")

if __name__ == '__main__':
    find_optimal_catalyst()