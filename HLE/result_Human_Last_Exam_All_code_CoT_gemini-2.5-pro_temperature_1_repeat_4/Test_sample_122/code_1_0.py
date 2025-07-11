def solve_tapestry_puzzle():
    """
    This function outlines the logical steps to solve the microbial metabolism puzzle
    and prints the final answer.
    """

    # Step 1: Analyze the final color.
    print("Step 1: The final color of the patch is orange.")
    print("Orange is a mixture of red and yellow pigments.")
    print("Therefore, the final patch must contain both red and yellow pigments to appear orange.")
    print("-" * 30)

    # Step 2: Determine the source of the red pigment and identify the mutation.
    print("Step 2: Trace the origin of the red pigment.")
    print("The metabolic pathway shows that red is an intermediate produced from yellow pigment by Enzyme A.")
    print("The reaction is: Yellow Pigment --(Enzyme A)--> Red Intermediate.")
    print("For red pigment to accumulate, the next step in the pathway must be blocked.")
    print("The next step is: Red Intermediate --(Enzyme B)--> Blue Intermediate.")
    print("To block this step, the microbe must be lacking Enzyme B.")
    print("Conclusion: The mutated enzyme is B.")
    print("-" * 30)

    # Step 3: Determine the original color of the patch.
    print("Step 3: Determine the origin of the yellow pigment in the final orange mixture.")
    print("The microbe is missing Enzyme B. The pathway that occurs is: Yellow Pigment --(Enzyme A)--> Red Intermediate.")
    print("If this reaction consumed all the yellow pigment, the final color would be red, not orange.")
    print("The presence of yellow pigment alongside the red pigment means that some of the original yellow pigment must remain.")
    print("Therefore, the patch must have been originally yellow.")
    print("-" * 30)

    # Step 4: Summarize the solution and print the final answer.
    final_enzyme = "B"
    original_colour = "yellow"
    print("Summary of the solution:")
    print(f"The original patch was {original_colour}.")
    print(f"The microbe is missing Enzyme {final_enzyme}.")
    print("Enzyme A converts some of the yellow pigment to red.")
    print(f"The process stops because Enzyme {final_enzyme} is missing, leaving a mixture of original yellow and new red pigment.")
    print("This mixture of yellow and red appears orange.")
    print("-" * 30)

    print("The final answer in the format <enzyme>-<colour> is:")
    print(f"{final_enzyme}-{original_colour}")

solve_tapestry_puzzle()
<<<B-yellow>>>