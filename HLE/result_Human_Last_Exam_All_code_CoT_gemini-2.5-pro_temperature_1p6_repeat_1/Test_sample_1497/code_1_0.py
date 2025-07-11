import sys

def solve_violin_puzzle():
    """
    This function determines the parameter groups corresponding to variations in playing a violin.

    The mapping is as follows:
    1. 'sul ponticello': Excites more high harmonics from the string. This changes the harmonic decay rate `ν`. Corresponds to group 'ii'.
    2. Bridge mute: Adds mass to the bridge, damping and lowering the frequency of body resonances, especially higher ones. This corresponds to group 'iv' (`μ, a₂, f₂`). Adding mass lowers frequency, so `f₂` goes 'down'.
    3. Helium room: Increases the speed of sound, which raises the frequency of the violin's internal air resonance (`f₁`). This corresponds to group 'iii'.
    4. E string vs A string: The fundamental frequency `F` of the open string is changed. This corresponds to group 'i'.

    The final answer is a comma-separated string of the results.
    """

    # Identify the group for each variation
    variation_1_group = "ii"
    variation_2_group = "iv"
    variation_3_group = "iii"
    variation_4_group = "i"

    # Identify the direction of change for the last member of the group for variation (2)
    variation_2_direction = "down"

    # Construct the final answer string
    final_answer = f"{variation_1_group},{variation_2_group},{variation_3_group},{variation_4_group},{variation_2_direction}"

    # Print the final answer
    print(final_answer)

solve_violin_puzzle()