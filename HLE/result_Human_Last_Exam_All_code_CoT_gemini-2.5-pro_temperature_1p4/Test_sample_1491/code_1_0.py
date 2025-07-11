def solve_replication_task():
    """
    Determines which AI capabilities are necessary for replicating the Visual Cliff and Swinging Room experiments.
    """
    # I. Goal-driven locomotion triggered by visual stimuli.
    # Necessary for the cliff (crawling to the goal) but not the room (standing still).
    answer_i = "cliff"

    # II. A functional analogue of cortical area MST in humans (optic flow processing).
    # Necessary for the room (processing moving wall illusion) and the cliff (processing motion parallax during self-motion).
    answer_ii = "both"

    # III. A perceptual mechanism for evaluating the relative size of elements in the visual scene (texture gradient).
    # Necessary for the cliff (a key monocular depth cue) but not the room.
    answer_iii = "cliff"

    # IV. A mechanism for resolving binocular disparity.
    # A powerful depth cue for the cliff, but not necessary for the room (optic flow works monocularly).
    answer_iv = "cliff"

    # V. A mechanism for controlling the effectors that maintain posture.
    # The primary system under investigation in the room experiment. While needed for crawling in the cliff,
    # it is the core measured component of the room paradigm.
    answer_v = "room"

    # Combine the answers in order, separated by hyphens.
    final_answer = "-".join([answer_i, answer_ii, answer_iii, answer_iv, answer_v])
    
    print(final_answer)

solve_replication_task()