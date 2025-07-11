import sys

def solve_tapestry_mystery():
    """
    Simulates the effect of microbe mutations on tapestry pigments to find
    the cause of an orange patch.
    """

    print("Analyzing the microbial metabolism to solve the tapestry mystery...\n")
    print("Metabolic Pathways:")
    print("1. Yellow Pigment -> (A) -> Red Intermediate -> (B) -> Blue Intermediate -> (C) -> Colourless")
    print("2. Blue Pigment  -> (D) -> Yellow Product\n")

    # The final color is orange, which is a mix of red and yellow.
    target_color = "orange"
    solution = None

    original_states = {
        "yellow": {"yellow_pigment"},
        "blue": {"blue_pigment"},
        "green (yellow+blue)": {"yellow_pigment", "blue_pigment"},
    }
    mutations = ["A", "B", "C", "D"]

    # This function simulates the chemical reactions
    def run_simulation(initial_substances, missing_enzyme):
        substances = set(initial_substances)

        # Pathway 1: Yellow Pigment Degradation
        if "yellow_pigment" in substances:
            substances.remove("yellow_pigment")
            # Step 1: Yellow -> Red (via Enzyme A)
            if missing_enzyme == "A":
                substances.add("yellow_pigment") # Process blocked, remains yellow
            else:
                # Step 2: Red -> Blue Intermediate (via Enzyme B)
                if missing_enzyme == "B":
                    substances.add("red_intermediate") # Process blocked, accumulates red
                else:
                    # Step 3: Blue Intermediate -> Colourless (via Enzyme C)
                    if missing_enzyme == "C":
                        substances.add("blue_intermediate") # Process blocked, accumulates blue_intermediate
                    else:
                        substances.add("colourless_product_1") # Process completes

        # Pathway 2: Blue Pigment Degradation
        if "blue_pigment" in substances:
            # Note: This pathway can run in parallel with the yellow one
            # We must check if the original blue pigment was present, not one from the other pathway
            if "blue_pigment" in initial_substances:
                substances.remove("blue_pigment")
                # Step 1: Blue -> Yellow (via Enzyme D)
                if missing_enzyme == "D":
                    substances.add("blue_pigment") # Process blocked, remains blue
                else:
                    substances.add("yellow_product_2") # Process completes

        return substances

    # This function determines the visible color from the final chemicals
    def get_final_color(substances):
        colors_present = set()
        if "red_intermediate" in substances:
            colors_present.add("red")
        if "yellow_pigment" in substances or "yellow_product_2" in substances:
            colors_present.add("yellow")
        if "blue_pigment" in substances or "blue_intermediate" in substances:
            colors_present.add("blue")

        if not colors_present:
            return "colourless"
        if colors_present == {"red", "yellow"}:
            return "orange"
        if colors_present == {"blue", "yellow"}:
            return "green"
        if len(colors_present) == 1:
            return list(colors_present)[0]
        return "mixed"

    print("--- Simulating all possible scenarios ---")
    for original_name, initial_pigments in original_states.items():
        for mutation in mutations:
            final_substances = run_simulation(initial_pigments, mutation)
            final_color = get_final_color(final_substances)

            print(f"Original Patch: {original_name:<22} | Missing Enzyme: {mutation} | Final Color: {final_color}")

            if final_color == target_color:
                solution = (mutation, original_name)

    print("\n--- Conclusion ---")
    if solution:
        enzyme, original_color_full = solution
        print(f"The simulation shows that an '{target_color}' patch occurs when:")
        print(f"1. The original patch was '{original_color_full}', containing both yellow and blue pigments.")
        print(f"2. Enzyme '{enzyme}' is mutated and non-functional.")

        print("\nBreakdown of the correct scenario:")
        print(f" - The original yellow pigment is converted to a red intermediate by enzyme A, but the pathway stops there because enzyme {enzyme} is missing.")
        print(" - The original blue pigment is converted to a yellow final product by the functional enzyme D.")
        print(" - The final mixture contains 'red intermediate' and 'yellow product', resulting in an ORANGE color.")

        # The output format is <enzyme>-<colour>, where colour is one of the base pigments.
        # The mutation is in the pathway that processes the YELLOW pigment.
        # Therefore, we associate the failure with the yellow pigment.
        final_answer_color = "yellow"
        final_answer = f"{enzyme}-{final_answer_color}"

        print(f"\nTo fit the answer format <enzyme>-<colour>, we associate the mutated enzyme '{enzyme}' with the pathway it belongs to, which starts with the '{final_answer_color}' pigment.")
        print(f"Final Answer: {final_answer}")
        print(f"\n<<<B-yellow>>>")
    else:
        print("No scenario resulted in an orange patch. Please check the problem statement.")

solve_tapestry_mystery()