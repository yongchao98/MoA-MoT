def solve_cathedral_puzzle():
    """
    Solves the pipe puzzle by finding the total number of pipes
    based on the given constraints.
    """
    in_tune_pipes = 200

    # From the poem, we deduce two key constraints:
    # 1. The total number of pipes (T) must be a multiple of lcm(3, 5) = 15.
    # 2. The number of out-of-tune pipes (U) must be a multiple of lcm(7, 4) = 28.
    # We also know that T = U + in_tune_pipes.

    # We need to find the smallest integer T > 200 that satisfies the constraints.
    # We can search for T by checking multiples of 15, starting from the first one > 200.
    total_pipes = 210
    while True:
        out_of_tune_pipes = total_pipes - in_tune_pipes
        
        # Check if the number of out-of-tune pipes is a multiple of 28.
        if out_of_tune_pipes % 28 == 0:
            # We found the correct numbers that fit all conditions.
            break
        
        # Move to the next multiple of 15.
        total_pipes += 15

    # The final question asks how many pipes the tuner must find, which is half of the lost ones.
    tuner_finds = out_of_tune_pipes / 2

    print(f"The total number of pipes in the cathedral is {total_pipes}.")
    print(f"The number of pipes still in tune is {in_tune_pipes}.")
    print(f"This means the number of out-of-tune (lost) pipes is {total_pipes} - {in_tune_pipes} = {out_of_tune_pipes}.")
    print("\nThe question is: How many pipes must the tuner find if half are realigned?")
    print("This is half of the out-of-tune pipes.")
    print("\nFinal Calculation:")
    print(f"{out_of_tune_pipes} / 2 = {int(tuner_finds)}")

solve_cathedral_puzzle()
<<<140>>>