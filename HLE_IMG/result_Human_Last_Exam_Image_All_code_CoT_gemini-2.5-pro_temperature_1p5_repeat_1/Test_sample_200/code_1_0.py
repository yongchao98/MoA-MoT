def solve_geology_question():
    """
    Analyzes the geological transition and identifies plausible causes.

    The transition from Rock A (dark, fine-grained, likely deep-water shale)
    to Rock B (lighter, likely shallower-water deposits) indicates a marine regression,
    which is a fall in relative sea level.

    We evaluate the options based on whether they cause a sea level fall:
    i. Onset of glaciation: Traps water in ice sheets, causing global sea level to fall. Plausible.
    ii. Coastal subsidence: Land sinks, causing relative sea level to rise (transgression). Not plausible.
    iii. Global warming: Melts ice and causes thermal expansion, raising global sea level. Not plausible.
    iv. Decreased Mid-Ocean Ridge Activity: Ocean basins increase in volume, causing global sea level to fall. Plausible.

    Therefore, the plausible explanations are i and iv.
    """
    plausible_events = ['i', 'iv']
    answer = ",".join(plausible_events)
    print(answer)

solve_geology_question()
print("<<<i,iv>>>")