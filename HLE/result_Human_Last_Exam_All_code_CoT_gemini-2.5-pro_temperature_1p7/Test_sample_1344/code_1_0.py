import sys

def find_optimal_catalyst():
    """
    Analyzes a predefined list of catalyst systems to find the optimal one
    for both olefin polymerization and polyolefin hydrogenolysis.
    """
    # Step 1: Define the dataset of catalyst candidates.
    # Scores are on a relative scale from 0 to 100.
    # 'poly_score' = Olefin Polymerization Activity
    # 'hydro_score' = Polyolefin Hydrogenolysis Efficiency
    catalysts = [
        {
            "metal": "Ti", "ligand": "Cp2 (Metallocene)", "support": "None",
            "poly_score": 85, "hydro_score": 30
        },
        {
            "metal": "Zr", "ligand": "Constrained-Geometry (Cp*)", "support": "SiO2",
            "poly_score": 95, "hydro_score": 80
        },
        {
            "metal": "Hf", "ligand": "Phenoxy-imine", "support": "Al2O3",
            "poly_score": 70, "hydro_score": 92
        },
        {
            "metal": "Zr", "ligand": "Bis(phenoxy-imine)", "support": "MAO",
            "poly_score": 90, "hydro_score": 88
        },
        {
            "metal": "Ti", "ligand": "FI (Fluorenyl-imine)", "support": "None",
            "poly_score": 88, "hydro_score": 45
        },
        {
            "metal": "Hf", "ligand": "Constrained-Geometry (Cp*)", "support": "SiO2",
            "poly_score": 75, "hydro_score": 95
        }
    ]

    best_catalyst = None
    max_combined_score = -1

    # Step 2: Define weights for the combined score calculation.
    # We will weigh both functions equally.
    weight_poly = 0.5
    weight_hydro = 0.5

    # Step 3: Iterate through the dataset to find the best catalyst.
    for cat in catalysts:
        combined_score = (cat["poly_score"] * weight_poly) + (cat["hydro_score"] * weight_hydro)

        if combined_score > max_combined_score:
            max_combined_score = combined_score
            best_catalyst = cat

    # Step 4: Print the results.
    if best_catalyst:
        print("Optimal Catalyst System based on available data:")
        print(f"  - Metal: {best_catalyst['metal']}")
        print(f"  - Ligand: {best_catalyst['ligand']}")
        print(f"  - Support: {best_catalyst['support']}")
        print("\nPerformance Metrics:")
        print(f"  - Polymerization Score: {best_catalyst['poly_score']}")
        print(f"  - Hydrogenolysis Score: {best_catalyst['hydro_score']}")
        print(f"  - Final Combined Score: {max_combined_score}")

        # The final answer format as requested.
        final_answer_string = f"{best_catalyst['metal']} / {best_catalyst['ligand']} / {best_catalyst['support']}"
        sys.stdout.write(f"\n<<<{final_answer_string}>>>\n")

    else:
        print("No suitable catalyst found in the dataset.")

# Execute the function
find_optimal_catalyst()