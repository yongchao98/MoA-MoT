def solve_psychology_replication_task():
    """
    Determines the necessary AI capabilities for replicating classic psychology experiments.

    This function analyzes five AI capabilities and decides whether each is necessary
    for replicating the "visual cliff" experiment, the "swinging room" experiment,
    "both", or "neither".

    The analysis is as follows:
    I. Goal-driven locomotion: Essential for the cliff paradigm where the subject must be motivated to move. Not for the room. -> cliff
    II. MST analogue (optic flow processing): Essential for the room paradigm which is entirely based on optic flow. Not essential for the cliff. -> room
    III. Relative size evaluation: Essential for the cliff as a key monocular depth cue (texture gradient). Not for the room. -> cliff
    IV. Binocular disparity: Essential for the cliff as a key binocular depth cue. Not for the room. -> cliff
    V. Postural control: Essential for both. The cliff requires it for locomotion, and the room requires it as the system being measured. -> both
    """
    
    # The answers for I through V in order
    capability_I = "cliff"
    capability_II = "room"
    capability_III = "cliff"
    capability_IV = "cliff"
    capability_V = "both"
    
    # Combine the answers into a single hyphen-separated string
    final_answer = f"{capability_I}-{capability_II}-{capability_III}-{capability_IV}-{capability_V}"
    
    print(final_answer)

solve_psychology_replication_task()