import sys
# Redirect print to a string buffer if not in an interactive environment
# This helps in environments where we want to capture the output as a whole.
# In a standard execution, this part has no effect.
try:
    __stdout__
except NameError:
    from io import StringIO
    sys.stdout = StringIO()


class Polysaccharide:
    """A simple model to describe a polysaccharide."""
    def __init__(self, monomer_name="Glucose"):
        self.type = "Polysaccharide"
        self.monomer = f"{monomer_name} (a monosaccharide)"
        self.backbone = "Chain of sugars linked by glycosidic bonds"
        self.elements = "Carbon (C), Hydrogen (H), Oxygen (O)"

    def describe(self):
        print(f"Polymer Type: {self.type}")
        print(f"Monomer Unit: {self.monomer}")
        print(f"Backbone Structure: {self.backbone}")
        print(f"Core Elements: {self.elements}")
        print("-" * 40)


class Polynucleotide:
    """A simple model to describe a polynucleotide."""
    def __init__(self, base="Adenine (A)"):
        self.type = "Polynucleotide"
        self.monomer = f"Nucleotide (sugar + PHOSPHATE + nitrogenous BASE)"
        self.backbone = "Alternating chain of sugars and PHOSPHATES"
        self.elements = "Carbon (C), Hydrogen (H), Oxygen (O), NITROGEN (N), PHOSPHORUS (P)"

    def describe(self):
        print(f"Polymer Type: {self.type}")
        print(f"Monomer Unit: {self.monomer}")
        print(f"Backbone Structure: {self.backbone}")
        print(f"Core Elements: {self.elements}")
        print("-" * 40)


def main():
    """Main function to run the comparison."""
    print("--- Comparing Polysaccharides and Polynucleotides ---\n")

    # Create and describe a polysaccharide (e.g., starch, which is a homopolysaccharide)
    starch = Polysaccharide(monomer_name="Glucose")
    starch.describe()

    # Create and describe a polynucleotide (e.g., poly-A DNA, a homopolynucleotide)
    poly_a_dna = Polynucleotide(base="Adenine (A)")
    poly_a_dna.describe()

    print("\nCONCLUSION:\n")
    print("Based on the structural differences, are polynucleotides polysaccharides?")
    final_answer = "No."
    print(f"The answer is a definitive '{final_answer}'")
    print("\nKey Differences:")
    print("1. Monomers: A polysaccharide's monomer is a simple sugar. A polynucleotide's monomer is a complex nucleotide, which contains a sugar, a phosphate group, and a nitrogenous base.")
    print("2. Backbone: A polysaccharide has a sugar-only backbone. A polynucleotide has a sugar-phosphate backbone.")
    print("3. Composition: Polysaccharides contain only C, H, and O. Polynucleotides also contain Nitrogen (in the bases) and Phosphorus (in the backbone).")

if __name__ == "__main__":
    main()

# Check if we buffered the output and print it now.
# In a standard execution, sys.stdout will be the original stdout object.
if hasattr(sys.stdout, 'getvalue'):
    output = sys.stdout.getvalue()
    sys.stdout = sys.__stdout__
    print(output)