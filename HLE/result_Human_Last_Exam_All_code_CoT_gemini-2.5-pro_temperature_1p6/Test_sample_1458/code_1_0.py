import sys

def solve_cathedral_riddle():
    """
    This function solves the word problem from "The Cathedral's Echo".
    It calculates the total number of pipes, the number of lost pipes,
    and finally, the number of pipes the tuner needs to find.
    """
    
    # The number of pipes that are still in tune.
    in_tune_pipes = 200
    
    # Let T be the total number of pipes.
    # The number of lost pipes (L) is T - 200.
    # The number of lost pipes is also described as (1/3)T + (2/5)T.
    # So, T - 200 = (1/3)T + (2/5)T
    # T - 200 = (5/15)T + (6/15)T
    # T - 200 = (11/15)T
    # T - (11/15)T = 200
    # (4/15)T = 200
    # T = 200 * 15 / 4
    
    # Calculate the total number of pipes (T)
    total_pipes = (200 * 15) / 4
    
    # Calculate the number of lost pipes (L)
    lost_pipes = total_pipes - in_tune_pipes
    
    # The final question asks how many pipes the tuner must find when half the lost ones are realigned.
    # This means the tuner needs to find the other half.
    tuner_finds = lost_pipes / 2

    # Print out the logic and the final answer in an equation format.
    print("Step 1: Find the total number of pipes (T).")
    print(f"The equation is: T - (1/3 * T + 2/5 * T) = {in_tune_pipes}")
    print(f"This simplifies to: (4/15) * T = {in_tune_pipes}")
    print(f"Solving for T: T = {in_tune_pipes} * 15 / 4 = {int(total_pipes)}")
    print(f"\nStep 2: Find the number of 'lost' pipes (L).")
    print(f"L = Total Pipes - In-Tune Pipes")
    print(f"L = {int(total_pipes)} - {in_tune_pipes} = {int(lost_pipes)}")
    print(f"\nStep 3: Calculate the number of pipes the tuner must find.")
    print(f"The tuner needs to find half of the lost pipes.")
    print(f"Final Equation: {int(lost_pipes)} / 2 = {int(tuner_finds)}")

solve_cathedral_riddle()
<<<275>>>