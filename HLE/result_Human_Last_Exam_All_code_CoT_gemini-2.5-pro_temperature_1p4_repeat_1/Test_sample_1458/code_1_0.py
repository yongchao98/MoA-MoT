def solve_cathedral_riddle():
    """
    Solves the Cathedral's Echo riddle by finding the smallest number of pipes
    that satisfies all the constraints given in the poem.
    """
    total_pipes = 0
    lost_pipes = 0

    # Start with a number that is a multiple of 15, as 'T' must be divisible by 3 and 5.
    # We can start checking from a reasonable lower bound.
    # T must be > 200 since there are 200 pipes in tune.
    # The smallest multiple of 15 > 200 is 210.
    current_t = 210

    while True:
        # The number of pipes "lost" (out of tune) is the total minus the pure ones.
        # This group L="lost their perfect pitch" must be divisible by 28 (lcm of 7 and 4).
        lost_pipes = current_t - 200
        if lost_pipes > 0 and lost_pipes % 28 == 0:
            total_pipes = current_t
            break
        # Increment to the next multiple of 15 and check again.
        current_t += 15

    # The final question asks for half of the "lost" pipes.
    tuner_finds = lost_pipes // 2

    print(f"The total number of pipes in the organ is: {total_pipes}")
    print(f"The number of pipes that still sing pure is: 200")
    print(f"The number of pipes that lost their perfect pitch is: {lost_pipes}")
    print("\nThe tuner must find the pipes to realign.")
    print("This is half of the lost pipes.")
    print(f"The final equation is: {lost_pipes} / 2 = {tuner_finds}")

solve_cathedral_riddle()