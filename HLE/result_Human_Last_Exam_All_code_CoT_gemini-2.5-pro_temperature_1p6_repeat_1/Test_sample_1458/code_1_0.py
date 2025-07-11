def solve_cathedral_echo():
    """
    Solves the Cathedral's Echo riddle by setting up and solving
    a system of equations based on the poem's text.
    """

    # We need to find a total number of pipes (T) that satisfies all conditions.
    # 1. T must be a multiple of 3 (for the 1/3 group).
    # 2. T must be a multiple of 5 (for the 2/5 group).
    # 3. The group that "lost their pitch" (1/3 * T) must be a multiple of 28
    #    (lcm of 7 and 4 for the subgroups).
    # This means T must be a multiple of lcm(5, 84) = 420.

    # We also know the total T = (in_tune) + (out_of_tune).
    # T = 200 + |union(L, M)| = 200 + |L| + |M| - |intersection(L, M)|
    # T = 200 + (1/3)T + (2/5)T - |intersection|
    # |intersection| = 200 + (11/15)T - T = 200 - (4/15)T
    # Since |intersection| >= 0, we have 200 >= (4/15)T, which means T <= 750.

    # The only multiple of 420 that is <= 750 is 420 itself.
    total_pipes = 420

    # The number of pipes "that lost their perfect pitch" is one-third of the total.
    pipes_lost = total_pipes / 3

    # The question asks how many pipes the tuner must find, which is half of the lost pipes.
    tuner_finds = pipes_lost / 2

    # Print the final equation as requested.
    print(f"The number of pipes that were lost is {int(pipes_lost)}.")
    print("The tuner must find half of these to realign them.")
    print(f"The final calculation is: {int(pipes_lost)} / 2 = {int(tuner_finds)}")

solve_cathedral_echo()
<<<70>>>