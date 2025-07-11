import sys

def solve_chemistry_puzzle():
    """
    This function provides the solution to the cubane synthesis puzzle.
    The four theoretical pairs of carbon atoms on the cubane product that can be substituted
    with a carboxylic acid are determined by the mechanism of the double Favorskii rearrangement.
    
    The reasoning leads to the following four pairs:
    1. The first COOH group is at C2, the second is at C4.
    2. The first COOH group is at C2, the second is at C8.
    3. The first COOH group is at C3, the second is at C4.
    4. The first COOH group is at C3, the second is at C8.
    
    This script will print these pairs in the requested format.
    """
    
    # The pairs of carbon atom numbers.
    # We sort the numbers within each pair for consistency.
    pair1 = (2, 4)
    pair2 = (2, 8)
    pair3 = (3, 4)
    pair4 = (3, 8)
    
    # We create a list of these pairs to format the output.
    all_pairs = [pair1, pair2, pair3, pair4]
    
    # Format the output string as per the user request, e.g., (a,b), (c,d), ...
    output_string = ", ".join([f"({p[0]},{p[1]})" for p in all_pairs])
    
    print(output_string)

solve_chemistry_puzzle()
# This is the final answer block. Do not edit.
sys.stdout.flush()