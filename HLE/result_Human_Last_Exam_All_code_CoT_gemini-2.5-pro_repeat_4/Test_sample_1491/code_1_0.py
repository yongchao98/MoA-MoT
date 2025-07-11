def solve_replication_puzzle():
    """
    This function determines which AI capabilities are necessary for replicating
    classic psychology experiments and prints the formatted answer.
    """
    # Analysis results for each capability from I to V
    capability_I = "both"
    capability_II = "room"
    capability_III = "cliff"
    capability_IV = "cliff"
    capability_V = "room"
    
    # Create a list of the answers in order
    answers = [
        capability_I,
        capability_II,
        capability_III,
        capability_IV,
        capability_V
    ]
    
    # Join the answers with a hyphen and print the final string
    final_answer_string = "-".join(answers)
    print(final_answer_string)

solve_replication_puzzle()