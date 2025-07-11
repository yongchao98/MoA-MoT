import math

def solve_cathedral_echo():
    """
    Solves the Cathedral's Echo riddle.
    """
    # Step 1 & 2: Define variables and knowns
    pipes_in_tune = 200

    # Step 3, 5, 6: Find the number of out-of-tune pipes (O)
    # We derived the equation 4*O + 15*I = 2200, where O is the number of
    # out-of-tune pipes and I is the intersection of the two out-of-tune groups.
    # O must be a multiple of lcm(4, 7) = 28.
    # We will loop through multiples of 28 to find a valid integer solution for I.
    
    out_of_tune_pipes = 0
    intersection_pipes = 0
    
    # We test multiples of 28 for O. Let's start from 28 and go up.
    k = 1
    while True:
        # O must be a multiple of 28
        O_candidate = 28 * k
        
        # From 4*O + 15*I = 2200, we get I = (2200 - 4*O) / 15
        # Check if the numerator is divisible by 15
        numerator = 2200 - 4 * O_candidate
        if numerator > 0 and numerator % 15 == 0:
            out_of_tune_pipes = O_candidate
            intersection_pipes = numerator // 15
            break # Found our solution
        k += 1
        if k > 100: # Safety break to prevent infinite loops
            print("Could not find a solution.")
            return

    # Step 7: Calculate all values for the final explanation
    total_pipes = pipes_in_tune + out_of_tune_pipes
    
    group_A_size = total_pipes // 3 # "One-third of pipes"
    group_B_size = (total_pipes * 2) // 5 # "two-fifths"
    
    octaves_rich = (out_of_tune_pipes * 3) // 7
    minor_scales = out_of_tune_pipes // 4
    discord_pipes = out_of_tune_pipes - octaves_rich - minor_scales
    
    # The final question: How many must the tuner find?
    # "When just half the lost realign"
    tuner_finds = out_of_tune_pipes // 2

    # Print the full breakdown of the solution
    print("The solution to the riddle is found through these steps:")
    print(f"1. The number of pipes in tune is {pipes_in_tune}.")
    print(f"2. The number of out-of-tune pipes ('O') must be a multiple of 28.")
    print(f"3. Solving the underlying equations, we find:")
    print(f"   - Total pipes in the cathedral: {total_pipes}")
    print(f"   - Pipes that 'fell out of tune' (1/3 of total): {group_A_size}")
    print(f"   - Pipes that 'caught the rising moon' (2/5 of total): {group_B_size}")
    print(f"   - Overlapping pipes in both groups: {intersection_pipes}")
    print("-" * 20)
    print(f"The total number of pipes that lost their perfect pitch is {out_of_tune_pipes}.")
    print("This is composed of:")
    print(f"   - {octaves_rich} pipes that found new octaves rich (3/7 of the lost).")
    print(f"   - {minor_scales} pipes that descended minor scales (1/4 of the lost).")
    print(f"   - {discord_pipes} pipes that were left in discord.")
    print("-" * 20)
    print("The final question is: How many must the tuner find when just half the lost realign?")
    print(f"The answer is half of the {out_of_tune_pipes} lost pipes.")
    print(f"Final Equation: {out_of_tune_pipes} / 2 = {tuner_finds}")

solve_cathedral_echo()
<<<140>>>