def solve_cathedral_riddle():
    """
    Solves the Cathedral's Echo riddle by finding the number of pipes based on divisibility constraints.
    """
    in_tune_pipes = 200
    out_of_tune_pipes = 0
    total_pipes = 0

    # We need to find the smallest number of out-of-tune pipes (O) that is a multiple of 28,
    # such that T = O + 200 is a multiple of 15.
    # We can test multiples of 28 until we find a valid solution.
    k = 1
    while True:
        # Out-of-tune pipes must be a multiple of 28 (lcm of 7 and 4)
        out_of_tune_pipes = 28 * k
        
        # Total pipes is the sum of in-tune and out-of-tune pipes
        total_pipes = in_tune_pipes + out_of_tune_pipes
        
        # Total pipes must be a multiple of 15 (lcm of 3 and 5)
        if total_pipes % 15 == 0:
            # We found the first set of numbers that satisfies all conditions.
            break
        
        k += 1

    # The final question is how many pipes need to be fixed, which is half the lost ones.
    pipes_to_fix = out_of_tune_pipes / 2

    print(f"Based on the riddle's constraints:")
    print(f"The number of pipes out of tune is a multiple of 28.")
    print(f"The total number of pipes is a multiple of 15.")
    print(f"The number of in-tune pipes is {in_tune_pipes}.")
    print("\nSolving for the smallest possible numbers:")
    print(f"Total pipes found (T): {total_pipes}")
    print(f"Pipes out of tune (O): {out_of_tune_pipes}")

    print("\nThe question asks how many pipes the tuner must find, which is half the lost (out-of-tune) pipes.")
    print("The final equation is:")
    print(f"{out_of_tune_pipes} / 2 = {int(pipes_to_fix)}")

solve_cathedral_riddle()
<<<140>>>