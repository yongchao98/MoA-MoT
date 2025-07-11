def solve_violin_puzzle():
    """
    This function determines the parameter group changes for different violin playing variations.

    The logic is as follows:
    1. 'sul ponticello' (bowing near the bridge) enhances high harmonics, which is controlled by the harmonic decay parameter 'ν' (group ii).
    2. A bridge mute adds mass, lowering the body's resonant frequencies, such as 'f_2' (group iv). The direction of change for frequency when adding mass is 'down'.
    3. Filling a room with helium increases the speed of sound, which raises the frequency of the violin's air-dependent resonances, such as 'f_1' (group iii).
    4. Playing on the E string instead of the A string changes the open string's fundamental frequency 'F' (group i).
    
    The final answer is formatted as a comma-separated string.
    """
    
    # Mapping of variations to parameter groups
    variation_1_group = "ii"  # sul ponticello -> ν
    variation_2_group = "iv"  # bridge mute -> f_2
    variation_3_group = "iii" # helium -> f_1
    variation_4_group = "i"   # E string -> F
    
    # Direction of change for the last member of the group for variation (2)
    # Adding mass (mute) lowers the resonant frequency f_2.
    direction_for_2 = "down"
    
    # Combine into the final answer format
    answer = f"{variation_1_group},{variation_2_group},{variation_3_group},{variation_4_group},{direction_for_2}"
    
    print(answer)

solve_violin_puzzle()