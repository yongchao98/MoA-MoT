def compare_macromolecules():
    """
    This script compares the basic building blocks of polysaccharides and polynucleotides
    to show why they are structurally different.
    """

    # Define the monomer for a polysaccharide
    polysaccharide_monomer = {
        "Name": "Monosaccharide (e.g., Glucose)",
        "Polymer Class": "Polysaccharide (e.g., Starch, Cellulose)",
        "Components": [
            "A simple sugar molecule containing only Carbon, Hydrogen, and Oxygen."
        ],
        "Bond Type": "Glycosidic Bond"
    }

    # Define the monomer for a polynucleotide
    polynucleotide_monomer = {
        "Name": "Nucleotide (e.g., Adenosine triphosphate)",
        "Polymer Class": "Polynucleotide (e.g., DNA, RNA)",
        "Components": [
            "A Five-Carbon Sugar (Deoxyribose or Ribose)",
            "A Phosphate Group",
            "A Nitrogenous Base"
        ],
        "Bond Type": "Phosphodiester Bond"
    }

    # Print the comparison
    print("--- Comparing the Building Blocks of Biopolymers ---\n")

    print("1. Polysaccharide Monomer:")
    for key, value in polysaccharide_monomer.items():
        if isinstance(value, list):
            print(f"   - {key}:")
            for item in value:
                print(f"     * {item}")
        else:
            print(f"   - {key}: {value}")

    print("\n------------------------------------------------------\n")

    print("2. Polynucleotide Monomer:")
    for key, value in polynucleotide_monomer.items():
        if isinstance(value, list):
            print(f"   - {key}:")
            for item in value:
                print(f"     * {item}")
        else:
            print(f"   - {key}: {value}")

    print("\n--- Conclusion ---")
    print("Although a polynucleotide's monomer contains a sugar, it also contains phosphate and a nitrogenous base, making it a nucleotide, not a monosaccharide.")
    print("Therefore, a polynucleotide is not a polysaccharide.")

if __name__ == '__main__':
    compare_macromolecules()