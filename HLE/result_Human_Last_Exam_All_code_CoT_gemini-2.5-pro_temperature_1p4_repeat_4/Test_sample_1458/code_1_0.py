import sys

def solve_cathedral_riddle():
    """
    This script solves the Cathedral's Echo riddle by determining the number of pipes a tuner must find.
    """
    # Step 1: Define the known values from the poem.
    in_tune_pipes = 200
    print(f"Based on the poem, the number of pipes that still sing pure is {in_tune_pipes}.")

    # Step 2: Establish the mathematical constraints.
    # The total number of pipes (T) must be a multiple of 3 and 5, so it must be a multiple of 15.
    # The number of lost pipes (L) must be a multiple of 7 and 4, so it must be a multiple of 28.
    # We also know that Total Pipes (T) = Lost Pipes (L) + In-Tune Pipes (I).
    # So, L + 200 must be a multiple of 15.
    print("The number of lost pipes, L, must be a multiple of 4 and 7 (i.e., of 28).")
    print(f"The total number of pipes, T, must be L + {in_tune_pipes}, and T must be a multiple of 3 and 5 (i.e., of 15).\n")
    print("Searching for the first number that satisfies these conditions...")

    lost_pipes = 0
    total_pipes = 0

    # Step 3: Find the number of lost pipes by checking multiples of 28.
    # We check multiples of 28 for L and see when L + 200 is a multiple of 15.
    # We can start our search from m=1.
    for m in range(1, 101): # A reasonable limit for the search
        potential_L = 28 * m
        potential_T = potential_L + in_tune_pipes
        
        if potential_T % 15 == 0:
            lost_pipes = potential_L
            total_pipes = potential_T
            print(f"Found a solution: If the number of lost pipes is {lost_pipes},")
            print(f"The total number of pipes is {lost_pipes} + {in_tune_pipes} = {total_pipes}.")
            print(f"Checking the condition: {total_pipes} is a multiple of 15 ({total_pipes} / 15 = {total_pipes//15}). This works!\n")
            break

    if lost_pipes == 0:
        print("Could not find a solution within the search limit.")
        sys.exit()

    # Step 4: Calculate the final answer based on the question.
    # The question is: "How many must the tuner find when just half the lost realign?"
    # This means we need to calculate half of the total lost pipes.
    pipes_to_find = lost_pipes // 2

    print("The riddle asks how many pipes the tuner must find if half of the lost pipes are fixed ('realigned').")
    print("This is equal to half of the total number of lost pipes.")
    print(f"The final calculation is: {lost_pipes} / 2 = {pipes_to_find}")


solve_cathedral_riddle()
<<<140>>>