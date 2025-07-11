import sys

def solve_tapestry_puzzle():
    """
    This script logically deduces the mutated enzyme and original tapestry colour.
    """
    
    # Step 1: Analyze the final observed colour.
    print("Step 1: The final colour of the patch is orange.")
    print("         Orange is a mixture of red and yellow colours.")
    print("-" * 30)

    # Step 2: Determine the source of the red colour.
    # Pathway 1: yellow pigment --(A)--> red intermediate --(B)--> blue intermediate --(C)--> colourless
    # Pathway 2: blue pigment --(D)--> yellow final product
    print("Step 2: Determine the source of the red component.")
    print("         According to the metabolic pathways, the only source of a red colour is the 'red intermediate'.")
    print("         This intermediate is produced from the original yellow pigment by Enzyme A.")
    print("         Therefore, the original patch must have contained yellow pigment, and Enzyme A must be functional.")
    print("-" * 30)

    # Step 3: Determine the source of the yellow colour.
    print("Step 3: Determine the source of the yellow component.")
    print("         Since the original colour was yellow, the yellow in the final orange mixture must be some of the original, unconverted yellow pigment.")
    print("-" * 30)

    # Step 4: Identify the specific mutation.
    print("Step 4: Identify the single enzyme mutation.")
    print("         For the 'red intermediate' to accumulate and not be converted further, the next step in its pathway must be blocked.")
    print("         The pathway is: ... -> red intermediate --(Enzyme B)--> blue intermediate -> ...")
    print("         A mutation in Enzyme B would block the pathway, causing the red intermediate to accumulate.")
    print("         The resulting mixture of original 'yellow pigment' and accumulated 'red intermediate' creates the colour orange.")
    print("-" * 30)

    # Step 5: Conclusion
    print("Conclusion:")
    print("The original colour of the patch was yellow.")
    print("The microbe is lacking Enzyme B.")
    
    mutated_enzyme = "B"
    original_colour = "yellow"
    
    print("\nFinal Answer:")
    print(f"{mutated_enzyme}-{original_colour}")

solve_tapestry_puzzle()

# The final answer in the required format.
sys.stdout.write("\n<<<B-yellow>>>\n")