import sys

def solve_violin_puzzle():
    """
    This function solves the violin parameter puzzle and prints the answer.
    
    Here is the reasoning:
    1.  Variation (1) 'sul ponticello': Bowing near the bridge excites more upper harmonics, creating a brighter sound.
        The term exp(-n/nu) controls the harmonic decay. A brighter sound means less decay, which corresponds to a larger nu.
        Thus, 'sul ponticello' primarily changes parameter group ii (nu).
    
    2.  Variation (2) 'with a bridge mute': A mute adds mass and damping to the bridge. This affects the resonances of the
        violin body, particularly those related to the bridge itself. The bridge's properties are critical for higher-frequency
        resonances. This corresponds to group iv (mu, a_2, f_2).
        
    3.  Variation (3) 'in a room filled with helium': The speed of sound in helium is ~2.8 times that in air. This
        dramatically increases the frequency of the main air resonance inside the violin body. This primary resonance is
        best represented by group iii (a_1, f_1).
        
    4.  Variation (4) 'on the E string': Playing on a different string (E instead of A) changes the fundamental
        frequency of the open string, F. This corresponds to group i (F).

    5.  Direction of change for (2):
        - Variation (2) is the mute, corresponding to group iv (mu, a_2, f_2).
        - The last member of the group is f_2, a resonance frequency.
        - Adding mass (the mute) to an oscillating system (the bridge) lowers its resonance frequency.
        - Therefore, the direction of change for f_2 is 'down'.
        
    The final sequence is: ii, iv, iii, i, down.
    """
    
    # Mapping each variation to its corresponding parameter group number
    # (using Roman numerals as requested)
    group_1 = "ii"  # sul ponticello -> nu
    group_2 = "iv"  # bridge mute -> mu, a2, f2
    group_3 = "iii" # helium -> a1, f1
    group_4 = "i"   # E string -> F
    
    # Direction of change for the last parameter of group 2 (f2)
    direction_2 = "down" # Adding mass lowers resonance frequency
    
    # Assemble the final answer string
    answer = f"{group_1},{group_2},{group_3},{group_4},{direction_2}"
    
    print(answer)

solve_violin_puzzle()