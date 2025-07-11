def solve_geology_question():
    """
    This function determines the plausible explanations for the geological transition shown in the image.

    The transition from the dark Rock A (likely a deep-water, anoxic deposit) to the lighter Rock B 
    (likely a shallower-water, more oxygenated deposit) indicates a marine regression, which is caused 
    by a fall in relative sea level.

    We evaluate the options:
    i. Onset of glaciation: Locks water in ice sheets, causing global sea-level FALL. -> Plausible.
    ii. Coastal subsidence: Land sinks, causing relative sea-level RISE. -> Not plausible.
    iii. Global warming: Melts ice and thermal expansion cause sea-level RISE. -> Not plausible.
    iv. Decreased Mid-Ocean Ridge Activity: Ocean basins expand in volume, causing global sea-level FALL. -> Plausible.
    
    The plausible options are 'i' and 'iv'.
    """
    
    # List of plausible explanations based on the analysis
    plausible_explanations = ['i', 'iv']
    
    # Format the answer as a comma-separated string with no spaces
    answer = ",".join(plausible_explanations)
    
    # Print the final answer
    print(answer)

solve_geology_question()