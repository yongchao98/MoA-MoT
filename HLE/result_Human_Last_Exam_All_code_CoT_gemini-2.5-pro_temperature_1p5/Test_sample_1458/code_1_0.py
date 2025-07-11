def solve_cathedral_echo():
    """
    Solves the riddle by finding numbers that satisfy all constraints.
    """
    pipes_in_tune = 200
    
    # We start k at 1 because there must be some pipes out of tune.
    k = 1
    while True:
        # Constraint 1: Out-of-tune pipes are a multiple of 28 (lcm of 7 and 4).
        pipes_out_of_tune = 28 * k
        
        # Total pipes are the sum of in-tune and out-of-tune pipes.
        total_pipes = pipes_out_of_tune + pipes_in_tune
        
        # Constraint 2: Total pipes must be a multiple of 15 (lcm of 3 and 5).
        if total_pipes % 15 == 0:
            # We found the first set of numbers that fits all constraints.
            
            # The final question asks for half the lost (out-of-tune) pipes.
            pipes_to_realign = pipes_out_of_tune / 2
            
            print(f"Based on the riddle's constraints:")
            print(f"The number of pipes still in tune is {pipes_in_tune}.")
            print(f"The number of pipes out of tune (O) must be a multiple of 28.")
            print(f"The total number of pipes (T) must be a multiple of 15.")
            print(f"Solving T = O + 200, we find:")
            print(f"Total Pipes = {total_pipes}")
            print(f"Pipes Out of Tune = {pipes_out_of_tune}")
            
            print("\nThe question is: How many must the tuner find when just half the lost realign?")
            print("The final equation is the number of out-of-tune pipes divided by two.")
            # Final output as requested
            print(f"{int(pipes_out_of_tune)} / 2 = {int(pipes_to_realign)}")
            
            # Return the final numerical answer
            return pipes_to_realign

        k += 1

# Execute the function to find the answer.
final_answer = solve_cathedral_echo()