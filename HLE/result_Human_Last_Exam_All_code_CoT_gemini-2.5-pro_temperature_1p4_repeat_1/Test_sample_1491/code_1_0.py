def solve_replication_paradigm():
    """
    This function determines which AI capabilities are necessary for replicating
    the Visual Cliff and Swinging Room experiments.
    """
    # The answers for capabilities I through V are determined based on the
    # core requirements of each psychological experiment.
    answers = [
        "both",  # I. Goal-driven locomotion triggered by visual stimuli.
        "room",  # II. A functional analogue of cortical area MST.
        "cliff", # III. A mechanism for evaluating relative size.
        "cliff", # IV. A mechanism for resolving binocular disparity.
        "room"   # V. A mechanism for controlling posture.
    ]

    # Join the answers with a hyphen as requested.
    final_answer = "-".join(answers)

    # Print the final formatted answer.
    print(final_answer)

solve_replication_paradigm()