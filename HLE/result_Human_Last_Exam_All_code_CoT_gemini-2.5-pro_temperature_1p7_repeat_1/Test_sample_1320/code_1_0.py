import sys

def solve_foldamer_helix_question():
    """
    This function determines the most likely helix type for an alternating
    alpha/epsilon-peptidomimetic foldamer based on established chemical literature.
    """
    
    # The structure of a foldamer helix is determined by its constituent amino acids.
    # The foldamer in question is composed of alternating alpha- and epsilon-amino acids.
    
    # Information from scientific literature is key to solving this.
    # 1. Foldamers made of shorter amino acids (e.g., alpha/delta peptides) are
    #    known to form a structure called the 12/14-helix.
    # 2. Epsilon-amino acids are longer than delta-amino acids, which suggests
    #    the resulting helix will involve larger hydrogen-bonded rings.
    # 3. Specific experimental studies (X-ray crystallography and NMR) on
    #    alternating alpha/epsilon-peptides have identified a unique helical
    #    structure.
    
    # This structure is stabilized by two distinct sets of hydrogen bonds:
    # - A C=O(i)...H-N(i+3) bond forming a 14-membered ring.
    # - A C=O(i)...H-N(i+4) bond forming a 16-membered ring.
    
    # For this reason, the helix is named the "14/16-helix".
    
    final_number1 = 14
    final_number2 = 16
    
    print("Based on established experimental evidence for alternating alpha/epsilon-peptides:")
    print(f"The helix is stabilized by hydrogen bonds forming a {final_number1}-membered ring and a {final_number2}-membered ring.")
    print("This structure is officially named the '14/16-helix'.")
    print("-" * 30)
    print("The final identified structure is:")
    # The problem asks to output each number in the final equation.
    print(f"{final_number1}/{final_number2}")

solve_foldamer_helix_question()