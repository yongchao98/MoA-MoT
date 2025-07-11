import math

def solve_cathedral_echo():
    """
    Solves the word problem from "The Cathedral's Echo".
    """
    # Known values from the riddle
    pipes_in_tune = 200
    fraction_out_of_tune_1 = 1/3
    fraction_out_of_tune_2 = 2/5
    
    # --- Step 1 & 2: Calculate the total fraction of pipes in and out of tune ---
    # The fractions are additive for the total set of pipes that fell out of tune.
    total_fraction_out_of_tune = fraction_out_of_tune_1 + fraction_out_of_tune_2
    fraction_in_tune = 1 - total_fraction_out_of_tune

    # --- Step 3: Calculate the total number of pipes ---
    # We know that the fraction of pipes in tune corresponds to 200 pipes.
    # So, (fraction_in_tune) * total_pipes = 200
    total_pipes = int(pipes_in_tune / fraction_in_tune)
    
    # --- Step 4: Calculate the total number of pipes out of tune ---
    pipes_out_of_tune = total_pipes - pipes_in_tune
    
    # --- Step 5: Find the final answer ---
    # The tuner must find half of the "lost" (out-of-tune) pipes.
    pipes_to_find = int(pipes_out_of_tune / 2)

    # --- Final Output ---
    # The fractions about octaves and scales are descriptive details not needed for the core calculation.
    
    print("Step 1: Determine the total number of pipes.")
    print(f"The fraction of pipes out of tune is 1/3 + 2/5 = 11/15.")
    print(f"The remaining fraction of pipes in tune is 1 - 11/15 = 4/15.")
    print(f"Given that 200 pipes are in tune, the total number of pipes is {pipes_in_tune} / (4/15) = {total_pipes}.")
    print("")
    print("Step 2: Determine the number of pipes the tuner must find.")
    print(f"The number of out-of-tune pipes ('the lost') is {total_pipes} (Total) - {pipes_in_tune} (In Tune) = {pipes_out_of_tune}.")
    print(f"The tuner must find half of these lost pipes.")
    print("")
    print("Final Equation:")
    print(f"{pipes_out_of_tune} / 2 = {pipes_to_find}")


solve_cathedral_echo()
<<<275>>>