import textwrap

def analyze_boronic_acid_synthesis():
    """
    Analyzes a problematic boronic acid synthesis and proposes a solution.
    """
    # --- Initial Problematic Reaction ---
    problem = {
        "starting_material": "2-bromo-4-chloro-1-iodobenzene",
        "reagents": {
            "n-BuLi": 1.05,  # equivalents
            "trimethyl borate": 5.0,  # equivalents
            "Solvent": "THF"
        },
        "temperature_C": -78,
        "observation": "Two different Boron (B) signals in NMR."
    }

    # --- Analysis ---
    print("--- Analysis of the Chemical Problem ---")
    print(f"Initial Observation: {problem['observation']}")
    print("\nExpected Product: (2-bromo-4-chlorophenyl)boronic acid should only have ONE boron signal.")
    
    explanation = """
    The presence of two boron signals is a known issue when a large excess of the borate ester is used.
    The primary product, the lithium boronate ate complex ([Ar-B(OMe)3]−Li+), can coordinate with a second molecule of the excess trimethyl borate.
    This creates a new complex with two chemically different boron environments, leading to two NMR signals.
    """
    print(textwrap.dedent(explanation))

    print("The key issue is the stoichiometry of the borating agent.")
    print(f"Initial Amount of Trimethyl Borate: {problem['reagents']['trimethyl borate']} equivalents (excessive).")

    # --- Solution ---
    print("\n--- Proposed Solution (Choice D) ---")
    solution_text = """
    To solve this, the excess of trimethyl borate should be reduced significantly.
    A smaller excess (e.g., 1.2 to 1.5 equivalents) is sufficient to trap the aryllithium without forming the problematic side-complex.
    """
    print(textwrap.dedent(solution_text))
    
    # --- Final Recommended Protocol ---
    corrected_reagents = problem['reagents'].copy()
    corrected_reagents['trimethyl borate'] = 1.5 # Corrected amount

    print("--- Final Recommended Reaction Equation ---")
    print("To obtain a single product, the corrected stoichiometry should be:")
    print(f"1.0 eq of {problem['starting_material']}")
    print(f"+ {corrected_reagents['n-BuLi']} eq of n-BuLi")
    print(f"+ {corrected_reagents['trimethyl borate']} eq of trimethyl borate")
    print(f"in {corrected_reagents['Solvent']} at {problem['temperature_C']} C")

if __name__ == "__main__":
    analyze_boronic_acid_synthesis()