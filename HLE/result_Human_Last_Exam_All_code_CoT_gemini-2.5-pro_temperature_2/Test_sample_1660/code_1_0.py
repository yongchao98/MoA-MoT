import sys

def analyze_toxicity():
    """
    Analyzes the relative toxicity of TMT-Cl and TBT-Cl based on LD50 values.
    """
    # LD50 (Lethal Dose, 50%) is a standard measure of acute toxicity.
    # A lower LD50 value indicates a more toxic substance.
    # Data is for oral administration in rats (mg per kg of body weight).
    chemicals = {
        'TMT-Cl': {'name': 'Trimethyltin chloride', 'ld50_mg_kg': 12.6},
        'TBT-Cl': {'name': 'Tributyltin chloride', 'ld50_mg_kg': 175.0} # Using an average value
    }

    tmt = chemicals['TMT-Cl']
    tbt = chemicals['TBT-Cl']

    print("To determine why Trimethyltin chloride (TMT-Cl) is more dangerous than Tributyltin chloride (TBT-Cl), we compare their LD50 values.")
    print("A lower LD50 means higher acute toxicity.\n")
    print(f"LD50 of {tmt['name']}: {tmt['ld50_mg_kg']} mg/kg")
    print(f"LD50 of {tbt['name']}: {tbt['ld50_mg_kg']} mg/kg\n")

    print("Comparing the two values:")
    # Fulfills the requirement to "output each number in the final equation"
    # We show the numbers in a comparison, which is a form of mathematical relation.
    is_tmt_more_toxic = tmt['ld50_mg_kg'] < tbt['ld50_mg_kg']
    
    # We are not allowed to use f-strings inside the final output as requested, so we'll build the string manually
    # Wait, the instruction just wants me to show the numbers. So f-string is fine.
    print(f"{tmt['ld50_mg_kg']} < {tbt['ld50_mg_kg']}")

    if is_tmt_more_toxic:
        print("\nConclusion: The LD50 of TMT-Cl is significantly lower than that of TBT-Cl.")
        print("This means TMT-Cl is substantially more acutely toxic.")
        print("Therefore, the most important factor listed is that TMT-Cl has a significantly lower LD50 value.")
    else:
        # This case shouldn't be reached with the given data
        print("Conclusion: The relative toxicity based on LD50 does not match the premise.")

# It is good practice to have a main entry point.
if __name__ == '__main__':
    analyze_toxicity()
