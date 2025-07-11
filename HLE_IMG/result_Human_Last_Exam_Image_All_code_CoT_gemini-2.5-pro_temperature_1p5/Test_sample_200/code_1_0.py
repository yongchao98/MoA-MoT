def solve_geology_question():
    """
    Analyzes the geological events that could cause the observed rock transition.

    The transition from Rock A (dark, fine-grained, deep-water deposit) to
    Rock B (lighter, shallower-water deposit) represents a marine regression,
    which is a relative fall in sea level.

    We evaluate the given options:
    i. Onset of glaciation: Causes a global sea-level fall. Plausible.
    ii. Coastal subsidence: Causes a relative sea-level rise. Not plausible.
    iii. Global warming: Causes a global sea-level rise. Not plausible.
    iv. Decreased Mid-Ocean Ridge Activity: Causes a global sea-level fall. Plausible.

    Therefore, the correct options are i and iv.
    """
    plausible_events = ['i', 'iv']
    answer = ",".join(plausible_events)
    print("The plausible explanations for the transition from Rock A to Rock B are events that cause a relative sea-level fall (regression).")
    print("Analyzing the options:")
    print("i. Onset of glaciation -> Causes sea-level fall -> Plausible")
    print("ii. Coastal subsidence -> Causes sea-level rise -> Not Plausible")
    print("iii. Global warming -> Causes sea-level rise -> Not Plausible")
    print("iv. Decreased Mid-Ocean Ridge Activity -> Causes sea-level fall -> Plausible")
    print("\nFinal Answer:")
    print(answer)

solve_geology_question()