def solve_geology_question():
    """
    Analyzes the geological events that could lead to the observed rock strata transition.

    The transition from Rock A (likely deep-water shale) to Rock B (likely shallower-water
    sandstone or limestone) indicates a sea-level fall, known as a regression.
    We need to identify which of the given options would cause a global sea-level fall.

    i. Onset of glaciation: Locks ocean water in ice sheets, causing sea-level fall. (Plausible)
    ii. Coastal subsidence: Local sinking of land, causing relative sea-level rise. (Not plausible)
    iii. Global warming: Melts ice and causes thermal expansion, leading to sea-level rise. (Not plausible)
    iv. Decreased Mid-Ocean Ridge Activity: Increases ocean basin volume, causing sea-level fall. (Plausible)

    Therefore, the correct options are i and iv.
    """
    correct_options = ['i', 'iv']
    answer = ",".join(correct_options)
    print("The transition from Rock A to Rock B represents a regression (fall in sea level).")
    print("The plausible causes for a regression among the choices are:")
    print("i. Onset of glaciation (locks water in ice, lowering sea level)")
    print("iv. Decreased Mid-Ocean Ridge Activity (increases ocean basin volume, lowering sea level)")
    print("\nThe final answer is:")
    print(answer)

solve_geology_question()