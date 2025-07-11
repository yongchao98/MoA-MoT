import sys

def solve_violin_puzzle():
    """
    This function solves the violin parameter puzzle.

    The reasoning is as follows:
    1.  (sul ponticello): This bowing technique emphasizes high harmonics. In the formula, the strength of harmonics is controlled by the decay parameter nu (ν). Thus, this corresponds to group ii.
    2.  (bridge mute): A mute adds mass to the bridge, damping vibrations and particularly suppressing high-frequency resonances of the violin body. This primarily affects the parameters for higher resonances, group iv (μ, a_2, f_2). Adding mass to an oscillator lowers its frequency, so f_2 goes down.
    3.  (helium-filled room): The resonant frequencies of the air inside the violin body (f_m) are proportional to the speed of sound in the gas. Helium has a much higher speed of sound than air, so the air resonance frequencies (especially the main one, f_1) will increase significantly. This corresponds to group iii (a_1, f_1).
    4.  (on the E string): The parameter F is the fundamental frequency of the open string. Changing from the A string to the E string is a direct change of F. This corresponds to group i.

    Combining these results gives the final answer.
    """
    
    # Mapping of variations to parameter groups
    # (1) sul ponticello -> ii
    # (2) bridge mute -> iv
    # (3) helium -> iii
    # (4) E string -> i
    # Direction for f_2 in (2) -> down
    
    answer = "ii,iv,iii,i,down"
    
    print(answer)

solve_violin_puzzle()