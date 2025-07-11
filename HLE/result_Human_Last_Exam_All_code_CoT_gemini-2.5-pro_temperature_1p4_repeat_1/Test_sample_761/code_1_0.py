class Monosaccharide:
    """A simple class to represent a monosaccharide monomer, like glucose."""
    def __init__(self, name="Glucose"):
        self.name = name

    def __str__(self):
        return f"{self.name} (a simple sugar)"

class Nucleotide:
    """A simple class to represent a nucleotide monomer."""
    def __init__(self, base_name="Adenine", sugar_name="Deoxyribose"):
        self.base = f"Nitrogenous Base ({base_name})"
        self.sugar = f"Pentose Sugar ({sugar_name})"
        self.phosphate = "Phosphate Group"

    def __str__(self):
        return f"Nucleotide (composed of {self.sugar}, {self.phosphate}, and {self.base})"

class Polysaccharide:
    """A class to model and describe a polysaccharide."""
    def __init__(self, monomer, name="Starch"):
        self.name = name
        self.monomer = monomer

    def describe_structure(self):
        print(f"--- Describing: {self.name} (a Polysaccharide) ---")
        print(f"1. Monomer (Repeating Unit): {self.monomer}")
        print("2. Backbone Structure: A chain of sugars linked directly to each other (Sugar-Sugar-Sugar...)")
        print("3. Linking Bond: Glycosidic Bond")
        print("4. Elemental Composition: Primarily Carbon, Hydrogen, Oxygen (C, H, O)")
        print("-" * 40)

class Homopolynucleotide:
    """A class to model and describe a homopolynucleotide."""
    def __init__(self, monomer_nucleotide, name="Poly(A) DNA"):
        self.name = name
        self.monomer = monomer_nucleotide

    def describe_structure(self):
        print(f"--- Describing: {self.name} (a Homopolynucleotide) ---")
        print(f"1. Monomer (Repeating Unit): {self.monomer}")
        print("2. Backbone Structure: A repeating chain of sugar and phosphate (Sugar-Phosphate-Sugar-Phosphate...)")
        print("3. Linking Bond: Phosphodiester Bond")
        print("4. Elemental Composition: Carbon, Hydrogen, Oxygen, Nitrogen, and Phosphorus (C, H, O, N, P)")
        print("-" * 40)

# --- Main execution to compare the two structures ---

# 1. Create a representative polysaccharide (like a part of starch)
glucose_monomer = Monosaccharide("Glucose")
starch_polymer = Polysaccharide(glucose_monomer)
starch_polymer.describe_structure()

# 2. Create a representative homopolynucleotide (like a poly-Adenine strand of DNA)
adenine_nucleotide = Nucleotide("Adenine", "Deoxyribose")
poly_a_polymer = Homopolynucleotide(adenine_nucleotide)
poly_a_polymer.describe_structure()

# 3. Print the final conclusion
print("\nConclusion:")
print("Polysaccharides and homopolynucleotides are structurally different.")
print("They are built from different monomers (simple sugars vs. nucleotides),")
print("are linked by different chemical bonds (glycosidic vs. phosphodiester),")
print("and have different backbone compositions (sugar-only vs. sugar-phosphate).")
