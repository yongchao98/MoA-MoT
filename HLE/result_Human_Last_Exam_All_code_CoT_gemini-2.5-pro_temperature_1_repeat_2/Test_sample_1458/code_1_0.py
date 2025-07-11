import math

def solve_cathedral_echo():
    """
    Solves the Cathedral's Echo riddle step-by-step.
    """
    # Step 1 & 2: Establish constraints to find the smallest possible total number of pipes (T).
    # The number of pipes that "lost their perfect pitch" (L) is T/3.
    # L must be divisible by 7 and 4, so it must be divisible by lcm(7, 4) = 28.
    # This means T/3 must be a multiple of 28, so T must be a multiple of 3 * 28 = 84.
    # The number of pipes that "caught the rising moon" (M) is 2T/5, so T must be divisible by 5.
    # Therefore, T must be a multiple of lcm(84, 5).
    lcm_val = math.lcm(84, 5)
    
    # We test multiples of lcm_val (420, 840, ...) to find T. Let's start with the smallest.
    T = lcm_val # Assume T = 420
    
    # Step 3: Verify T=420 with the information from the poem.
    # T = I + L + M - |L intersect M|
    # T = 200 + T/3 + 2T/5 - |L intersect M|
    # T - T/3 - 2T/5 = 200 - |L intersect M|
    # T * (1 - 1/3 - 2/5) = 200 - |L intersect M|
    # T * (4/15) = 200 - |L intersect M|
    # intersection = 200 - T * (4/15)
    
    pipes_I = 200 # In tune pipes
    pipes_L = T / 3
    pipes_M = T * 2 / 5
    
    # Calculate the size of the intersection. A valid solution requires intersection >= 0.
    intersection = pipes_I - (T * (1 - 1/3 - 2/5))
    
    # Total pipes are the sum of in-tune pipes and the union of the two out-of-tune groups
    total_out_of_tune = pipes_L + pipes_M - intersection
    
    print(f"Based on the poem's constraints, the total number of pipes is {int(T)}.")
    print(f"The number of pipes that 'lost their perfect pitch' is {int(pipes_L)}.")
    print(f"The total number of out-of-tune pipes is {int(total_out_of_tune)}.")
    
    # Step 4 & 5: Answer the final question.
    # "How many must the tuner find / When just half the lost realign?"
    # "The lost" refers to group L, the pipes that lost their pitch.
    realigned_pipes = pipes_L / 2
    
    # The number to find is the remaining out-of-tune pipes.
    pipes_to_find = total_out_of_tune - realigned_pipes
    
    print("\nThe final question asks how many pipes the tuner must find.")
    print("This is the total number of out-of-tune pipes minus the ones that were realigned.")
    print("\nThe final calculation is:")
    print(f"{int(total_out_of_tune)} - {int(realigned_pipes)} = {int(pipes_to_find)}")

solve_cathedral_echo()
<<<150>>>