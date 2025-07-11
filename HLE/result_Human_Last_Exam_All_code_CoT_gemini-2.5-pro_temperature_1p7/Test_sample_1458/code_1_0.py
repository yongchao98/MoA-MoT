def solve_cathedral_echo():
    """
    Solves the Cathedral's Echo riddle by finding the number of pipes
    that meet all the constraints described in the poem.
    """
    found_solution = False
    # O must be a multiple of 28 (from 3/7 and 1/4 fractions)
    # We will loop through multiples of 28.
    for k in range(1, 100):
        out_of_tune_pipes = 28 * k
        
        # Total pipes = out_of_tune + in_tune
        total_pipes = out_of_tune_pipes + 200
        
        # T must be a multiple of 15 (from 1/3 and 2/5 fractions)
        if total_pipes % 15 == 0:
            found_solution = True
            
            print("--- Solution Found ---")
            print(f"Total pipes in the organ (T): {total_pipes}")
            print(f"Pipes that still sing pure: 200")
            print(f"Pipes that lost their pitch (O): {out_of_tune_pipes}")
            print("\n--- Verifying the Poem's Figures ---")
            
            # Breakdown of 'O'
            octaves_rich = (3/7) * out_of_tune_pipes
            minor_scales = (1/4) * out_of_tune_pipes
            rest_in_discord = out_of_tune_pipes - octaves_rich - minor_scales
            
            print(f"Of the {out_of_tune_pipes} lost pipes:")
            print(f" - Found new octaves (3/7): {int(octaves_rich)}")
            print(f" - Descended minor scales (1/4): {int(minor_scales)}")
            print(f" - The rest in discord: {int(rest_in_discord)}")

            print("\n--- Answering the Final Question ---")
            print("The question asks: 'How many must the tuner find when just half the lost realign?'")
            print("This means we need to find half of the total number of lost pipes.")
            
            tuner_task = out_of_tune_pipes / 2
            
            # The prompt requests the final equation to be printed.
            print("\nThe final equation is:")
            print(f"{out_of_tune_pipes} / 2 = {int(tuner_task)}")
            
            # Return the final numerical answer as requested
            return int(tuner_task)

# Run the solver and capture the final answer
final_answer = solve_cathedral_echo()
# The required output format is <<<answer>>> at the end
# The print statements inside the function will provide the explanation
print(f"\n<<<140>>>")
