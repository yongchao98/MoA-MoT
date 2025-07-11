import sys

def solve_violin_puzzle():
    """
    This function analyzes the effect of four variations on a violin's sound
    and maps each variation to a specific parameter group from a given model.
    It then assembles and prints the final answer in the required format.
    """

    # --- Reasoning for each choice ---

    # Variation 1: 'sul ponticello'
    # Playing near the bridge enhances high harmonics. This is controlled by the
    # harmonic decay parameter ν. A brighter sound means larger ν.
    # This is a change in group ii.
    variation_1_group = "ii"

    # Variation 2: 'with a bridge mute'
    # A mute adds mass to the bridge, damping high-frequency resonances.
    # This affects the body resonance parameters, particularly the higher ones (group iv).
    # Adding mass lowers the resonant frequency (f₂).
    variation_2_group = "iv"
    direction_for_2 = "down"

    # Variation 3: 'in a room filled with helium'
    # Air resonance frequencies depend on the speed of sound, which is much higher
    # in helium. This shifts air resonances (like the one modeled by f₁) upwards.
    # This is a change in group iii.
    variation_3_group = "iii"

    # Variation 4: 'on the E string' (instead of A string)
    # The parameter F is the fundamental frequency of the open string.
    # Switching from the A string to the E string directly changes F.
    # This is a change in group i.
    variation_4_group = "i"

    # --- Assemble and print the final answer ---

    # The final answer is a comma-separated string of the identified groups
    # and the direction for the second variation.
    final_answer = ",".join([
        variation_1_group,
        variation_2_group,
        variation_3_group,
        variation_4_group,
        direction_for_2
    ])

    print(final_answer)

# Execute the function to solve the puzzle.
solve_violin_puzzle()