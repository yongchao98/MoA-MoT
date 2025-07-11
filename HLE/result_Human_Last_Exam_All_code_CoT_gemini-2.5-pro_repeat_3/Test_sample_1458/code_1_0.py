def solve_cathedral_echo():
    """
    Solves the Cathedral's Echo riddle by setting up and solving a system of constraints.
    """
    # From the analysis, the smallest possible total number of pipes (T) that satisfies
    # all fractional conditions (T divisible by 3 and 5; T/3 divisible by 7 and 4)
    # is the least common multiple of 15 and 84, which is 420.
    total_pipes = 420

    # The number of pipes given as still in tune.
    in_tune_pipes = 200

    # Calculate the size of the group that "lost their perfect pitch" (L).
    lost_pitch_pipes = total_pipes // 3
    
    # The question asks how many pipes the tuner must find. This is identified as
    # half of the "lost" pipes, which is the group L.
    pipes_to_find = lost_pitch_pipes // 2

    print(f"First, we determined the total number of pipes in the cathedral.")
    print(f"Based on the fractions in the poem, the smallest possible total is {total_pipes} pipes.")
    print(f"The group of pipes described as having 'lost their perfect pitch' is one-third of this total.")
    print(f"Number of 'lost' pipes = {total_pipes} / 3 = {lost_pitch_pipes}")
    print(f"The question is how many of these {lost_pitch_pipes} pipes the tuner must find if half are realigned.")
    print(f"The tuner must find the other half, which are still out of tune.")
    print(f"The final equation is:")
    print(f"{lost_pitch_pipes} / 2 = {pipes_to_find}")

solve_cathedral_echo()
<<<70>>>