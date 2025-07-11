import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

def identify_foldamer_helix():
    """
    Identifies the likely helical structure of a given peptidomimetic foldamer
    based on established chemical principles.
    """

    # A knowledge base mapping foldamer compositions to their known helical structures.
    # The key is a tuple: (monomer1_type, monomer2_type, arrangement)
    foldamer_helix_database = {
        ("alpha", "epsilon", "alternating"): {
            "name": "14/16-helix",
            "ring_size_1": 14,
            "ring_size_2": 16,
            "choice": "H",
            "explanation": (
                "Foldamers composed of alternating alpha-amino acids and epsilon-amino acids "
                "are known to form a stable, predictable structure called a 14/16-helix. "
                "This helix is stabilized by two different, co-existing types of hydrogen bonds."
            )
        }
    }

    # Define the peptidomimetic from the user's query
    peptide_composition = ("alpha", "epsilon", "alternating")

    # Look up the helix type in our knowledge base
    if peptide_composition in foldamer_helix_database:
        helix_info = foldamer_helix_database[peptide_composition]
        
        print(f"The foldamer is an alternating copolymer of alpha- and epsilon-amino acids.")
        print(f"Based on established research in foldamer chemistry, this composition forms a specific structure: the {helix_info['name']}.")
        print("\nExplanation of the helix structure:")
        print(helix_info['explanation'])
        
        # The final answer format requires printing each number in the final equation/name.
        ring1 = helix_info['ring_size_1']
        ring2 = helix_info['ring_size_2']
        
        print(f"\n- The first hydrogen bond pattern forms a ring containing {ring1} atoms.")
        print(f"- The second hydrogen bond pattern forms a ring containing {ring2} atoms.")
        print(f"\nTherefore, the correct notation for this helix is: {ring1}/{ring2}")

    else:
        print("The specific foldamer composition is not found in the knowledge base.")

# Execute the function to generate the explanation
identify_foldamer_helix()

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output_str = captured_output.getvalue()

# Print the captured output to the actual console
print(output_str)

# Final answer selection based on the analysis
final_answer = "H"
print(f'<<<H>>>')