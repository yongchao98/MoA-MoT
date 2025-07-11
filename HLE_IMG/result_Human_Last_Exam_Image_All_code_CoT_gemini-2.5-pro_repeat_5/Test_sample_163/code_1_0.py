import sys
from io import StringIO

# A simple class to represent chemical compounds for clarity
class Compound:
    def __init__(self, name, smiles):
        self.name = name
        self.smiles = smiles

    def __str__(self):
        return f"Name: {self.name}\nSMILES: {self.smiles}"

def identify_reaction_products():
    """
    Identifies and prints the structures of products A and B from the reaction
    of styrene with tert-butyl peroxybenzoate catalyzed by Fe(OTf)3.
    """
    # Explanation of the reaction mechanism leading to two products
    explanation = """
The reaction involves the generation of a tert-butoxy radical, which has two competing fates:

1.  Pathway to Product A: The tert-butoxy radical adds to styrene, and the resulting benzylic radical is trapped by a benzoate group.
2.  Pathway to Product B: The tert-butoxy radical decomposes into a methyl radical and acetone. The methyl radical then adds to styrene, and the resulting benzylic radical is trapped by a benzoate group.
"""

    # Define the two major products
    product_A = Compound(
        name="2-(tert-butoxy)-1-phenylethyl benzoate",
        smiles="CC(C)(C)OCC(c1ccccc1)OC(=O)c2ccccc2"
    )

    product_B = Compound(
        name="1-phenylpropyl benzoate",
        smiles="CCC(c1ccccc1)OC(=O)c2ccccc2"
    )
    
    # Capture original stdout
    original_stdout = sys.stdout
    # Create a string buffer
    string_buffer = StringIO()
    # Redirect stdout to the buffer
    sys.stdout = string_buffer

    # Print the results to the buffer
    print("Based on the reaction mechanism, the two major products A and B are:\n")

    print("--- Product A ---")
    print(product_A)
    print("\n--- Product B ---")
    print(product_B)
    
    # Get the content of the buffer
    output_string = string_buffer.getvalue()
    # Restore original stdout
    sys.stdout = original_stdout
    # Print the captured output to the console
    print(output_string)


if __name__ == "__main__":
    identify_reaction_products()