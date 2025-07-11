import sys

def solve_mating_age_puzzle():
    """
    This function explains and prints the estimated mating age for two insect species
    based on their life cycles and the definition of 'eclosion'.
    """

    # --- Species 1: Tridactylophagus tartari (Wasp) ---
    species1_name = "Tridactylophagus tartari"
    species1_age_value = 3
    species1_age_unit = "days"
    species1_lifecycle = "Holometabolous (complete metamorphosis)"
    species1_eclosion_meaning = "Emergence of the adult from the pupa."

    # --- Species 2: Periplaneta americana (Cockroach) ---
    species2_name = "Periplaneta americana"
    species2_age_value = 6
    species2_age_unit = "months"
    species2_lifecycle = "Hemimetabolous (incomplete metamorphosis)"
    species2_eclosion_meaning = "Hatching from the egg."

    # Print the step-by-step reasoning
    print(f"Analysis for {species1_name}:")
    print(f"  - Life Cycle: {species1_lifecycle}")
    print(f"  - Eclosion Definition: {species1_eclosion_meaning}")
    print("  - Reasoning: Adult wasps mature very quickly after emerging to find a mate. A few days is a standard estimate for this to occur.")
    print("-" * 40)
    
    print(f"Analysis for {species2_name}:")
    print(f"  - Life Cycle: {species2_lifecycle}")
    print(f"  - Eclosion Definition: {species2_eclosion_meaning}")
    print("  - Reasoning: The age since hatching includes the long nymphal development period, which averages 6-12 months before the cockroach becomes a sexually mature adult.")
    print("-" * 40)
    
    # Print the final conclusion
    print("Final Estimate:")
    print(f"The average age for a male {species1_name} is {species1_age_value} {species1_age_unit}.")
    print(f"The average age for a male {species2_name} is {species2_age_value} {species2_age_unit}.")

solve_mating_age_puzzle()
<<<C>>>