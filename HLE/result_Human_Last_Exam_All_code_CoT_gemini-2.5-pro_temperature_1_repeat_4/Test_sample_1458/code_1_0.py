def solve_cathedral_riddle():
    """
    Solves the riddle by finding numbers that satisfy the given constraints.
    
    - The number of out-of-tune pipes (O) must be a multiple of 28 (from lcm(7,4)).
    - The total number of pipes (T) must be a multiple of 15 (from lcm(3,5)).
    - The number of in-tune pipes is 200, so T = O + 200.
    """
    
    # Start with the first possible number of out-of-tune pipes.
    out_of_tune_pipes = 28
    
    while True:
        # Calculate the corresponding total number of pipes.
        total_pipes = out_of_tune_pipes + 200
        
        # Check if the total number of pipes is divisible by 15.
        if total_pipes % 15 == 0:
            # If it is, we have found our solution.
            break
        
        # If not, move to the next multiple of 28.
        out_of_tune_pipes += 28
        
    # The question asks for half the number of lost (out-of-tune) pipes.
    tuner_must_find = out_of_tune_pipes / 2
    
    print(f"Based on the poem's constraints:")
    print(f"The number of out-of-tune ('lost') pipes is: {out_of_tune_pipes}")
    print(f"The total number of pipes in the cathedral is: {total_pipes}")
    print("\nThe question asks how many pipes the tuner must find, which is half of the lost pipes.")
    print(f"The final calculation is: {out_of_tune_pipes} / 2 = {int(tuner_must_find)}")

solve_cathedral_riddle()