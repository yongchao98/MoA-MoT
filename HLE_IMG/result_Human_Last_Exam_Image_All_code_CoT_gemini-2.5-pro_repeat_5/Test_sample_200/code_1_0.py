def solve_geology_problem():
    """
    Analyzes the geological sequence and determines the plausible causal events.

    The image shows a transition from a lower rock layer (A) to an upper rock layer (B).
    - Rock A: Dark, finely-layered, suggesting a deep-water, low-energy environment (e.g., offshore shale).
    - Rock B: Lighter, more massive, suggesting a shallower-water, higher-energy environment (e.g., nearshore sandstone or limestone).

    This transition from deep to shallow water deposits is called a regression, which is caused by a fall in relative sea level. We must identify which events cause a sea-level fall.

    Let's evaluate the options:
    i. Onset of glaciation: Locks ocean water in ice sheets, causing global sea level to fall. This causes a regression. (Plausible)
    ii. Coastal subsidence: The land sinks, causing a relative sea-level rise (transgression). (Not plausible)
    iii. Global warming: Melts ice sheets, causing global sea level to rise (transgression). (Not plausible)
    iv. Decreased Mid-Ocean Ridge Activity: Ocean basins increase in volume, causing global sea level to fall. This causes a regression. (Plausible)

    Therefore, the plausible explanations are i and iv.
    """
    
    # List of all options
    options = {
        'i': "Onset of glaciation",
        'ii': "Coastal subsidence",
        'iii': "Global warming",
        'iv': "Decreased Mid-Ocean Ridge Activity"
    }

    # Identifying the correct options based on geological reasoning
    correct_options = ['i', 'iv']

    # Formatting the answer as a comma-separated string with no spaces
    answer = ",".join(correct_options)
    
    print(answer)

solve_geology_problem()
<<<i,iv>>>