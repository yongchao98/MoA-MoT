def solve_cathedral_puzzle():
    """
    Solves the word puzzle about the cathedral pipes.
    """
    
    out_of_tune_pipes = 0
    total_pipes = 0
    
    # We are looking for the number of out-of-tune pipes, let's call it 'O'.
    # From the poem, 'O' must be a multiple of lcm(7, 4), which is 28.
    # The total number of pipes, 'T', must be a multiple of lcm(3, 5), which is 15.
    # We also know T = O + 200 (the number of pipes that still sing pure).
    
    # We'll test multiples of 28 for 'O' until we find one where O + 200 is a multiple of 15.
    o_candidate = 28
    while True:
        t_candidate = o_candidate + 200
        if t_candidate % 15 == 0:
            out_of_tune_pipes = o_candidate
            total_pipes = t_candidate
            break
        o_candidate += 28

    # The final question asks for half of the lost pipes.
    pipes_to_find = out_of_tune_pipes / 2
    
    print("Based on the poem's constraints:")
    print(f"The total number of pipes must be a multiple of 15.")
    print(f"The number of out-of-tune pipes ('the lost') must be a multiple of 28.")
    print(f"There are 200 pipes that are still in tune.")
    print(f"The first number to satisfy these conditions is {out_of_tune_pipes} out-of-tune pipes, for a total of {total_pipes} pipes.")
    print("\nThe question is: How many pipes must the tuner find if half the lost are realigned?")
    print("This requires calculating half of the out-of-tune pipes.")
    print("\nThe final equation is:")
    print(f"{out_of_tune_pipes} / 2 = {int(pipes_to_find)}")

solve_cathedral_puzzle()