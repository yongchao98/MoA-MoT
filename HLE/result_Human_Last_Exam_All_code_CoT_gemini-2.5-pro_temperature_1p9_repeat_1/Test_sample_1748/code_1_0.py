def solve_rope_riddle():
    """
    This function solves the magical rope riddle step by step.
    """
    initial_length = 60
    shrink_factor = 0.75 # 100% - 25%

    print("Step 1: The First Cut")
    print(f"Start with the initial rope of {initial_length}cm.")
    
    # We propose a cut that leads to the solution.
    cut1_long = 40
    cut1_short = 20
    print(f"Cut the rope into two pieces: {cut1_long}cm and {cut1_short}cm.")
    print(f"Check: {cut1_long} + {cut1_short} = {cut1_long + cut1_short}cm.")
    print("-" * 20)
    
    print("Step 2: The First Shrink")
    print("The longer portion from the cut (40cm) shrinks by 25%.")
    shrunk1_length = int(cut1_long * shrink_factor)
    print(f"Calculation: {cut1_long} * (1 - 0.25) = {shrunk1_length}")
    print(f"After shrinking, the two pieces are now {shrunk1_length}cm and {cut1_short}cm.")
    print("-" * 20)

    # Determine the longer portion after the first cut
    longer_portion_after_cut1 = max(shrunk1_length, cut1_short)
    
    print("Answering the question:")
    print("What will be the length of the longer portion after the first cut?")
    print(f"The two pieces are {shrunk1_length}cm and {cut1_short}cm. The longer one is {longer_portion_after_cut1}cm.")
    print("-" * 20)

    print("Step 3: The Second Cut (Proof of a valid solution)")
    # We take the longer piece from the first cut for our second cut
    rope_for_cut2 = longer_portion_after_cut1
    print(f"To get a 15cm piece, we take the {rope_for_cut2}cm rope for the second cut.")
    # To get 15 after shrinking, the original piece must be 15 / 0.75 = 20
    cut2_long = 20
    cut2_short = rope_for_cut2 - cut2_long
    print(f"Cut the {rope_for_cut2}cm rope into {cut2_long}cm and {cut2_short}cm pieces.")
    print("-" * 20)

    print("Step 4: The Second Shrink")
    print("The longer portion from this cut (20cm) shrinks by 25%.")
    shrunk2_length = int(cut2_long * shrink_factor)
    print(f"Calculation: {cut2_long} * (1 - 0.25) = {shrunk2_length}")
    print(f"The final pieces are {shrunk2_length}cm and {cut2_short}cm.")
    print(f"We have successfully obtained a {shrunk2_length}cm piece of rope.")

solve_rope_riddle()
<<<30>>>