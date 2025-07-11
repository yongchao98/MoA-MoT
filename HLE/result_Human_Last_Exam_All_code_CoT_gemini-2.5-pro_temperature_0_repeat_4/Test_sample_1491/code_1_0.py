def solve_paradigm_puzzle():
    """
    This function determines the required capabilities for replicating two classic psychology experiments.
    I. Goal-driven locomotion triggered by visual stimuli. -> both
    II. A functional analogue of cortical area MST in humans. -> both
    III. A perceptual mechanism for evaluating the relative size of elements in the visual scene. -> cliff
    IV. A mechanism for resolving binocular disparity. -> cliff
    V. A mechanism for controlling the effectors that maintain posture. -> both
    """
    
    # The answers for I through V are determined based on the analysis of each experiment.
    answers = ["both", "both", "cliff", "cliff", "both"]
    
    # The final output is formatted as a hyphen-separated string.
    final_answer = "-".join(answers)
    
    print(final_answer)

solve_paradigm_puzzle()