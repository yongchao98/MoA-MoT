import math

def solve_cathedral_riddle():
    """
    This script solves the word puzzle from the poem "The Cathedral's Echo".
    It calculates the total number of pipes and the number of pipes a tuner must find.
    """
    
    # Step 1: Define known variables from the poem.
    pure_pipes = 200
    print(f"Step 1: Identify the knowns.")
    print(f"There are {pure_pipes} pipes that still sing pure.")
    print("Let L be the number of pipes that 'lost their perfect pitch'.")
    print("Let T be the total number of pipes.")
    print("The total number of pipes is the sum of pure and lost pipes: T = L + 200.\n")

    # Step 2: Determine the properties of L, the number of lost pipes.
    print("Step 2: Analyze the composition of the lost pipes.")
    print("Of the lost pipes (L), 3/7 went to new octaves and 1/4 descended.")
    # For these fractions to result in a whole number of pipes, L must be a multiple of 7 and 4.
    lcm_denom = math.lcm(7, 4)
    print(f"Therefore, L must be a multiple of the Least Common Multiple of 7 and 4, which is {lcm_denom}.\n")
    
    # Step 3: Use the origin of the dissonance to find the final relationship.
    print("Step 3: Establish the equation linking all variables.")
    print("The lost pipes (L) are the union of two groups: T/3 (from thunder) and 2T/5 (from the moon).")
    print("The relationship is: L = (T/3) + (2T/5) - Intersection.")
    print("Substituting T = L + 200 and rearranging gives the equation: 2200 - 4*L = 15 * Intersection.")
    print("This means (2200 - 4*L) must be a non-negative number divisible by 15.\n")
    
    # Step 4: Solve for L by iterating through its possible values (multiples of 28).
    print("Step 4: Find the correct value for L by testing multiples of 28.")
    lost_pipes = 0
    # We can iterate through k, where L = 28 * k
    for k in range(1, 20): # k must be < 19.64 from 2200-4*28k > 0
        l_candidate = lcm_denom * k
        # Check if the condition (2200 - 4L) is a multiple of 15 is met.
        if (2200 - 4 * l_candidate) % 15 == 0:
            lost_pipes = l_candidate
            print(f"The first multiple of {lcm_denom} that fits the equation is {lost_pipes}.")
            break
    
    print(f"So, the number of lost pipes is {lost_pipes}.\n")

    # Step 5: Verify the total number of pipes.
    total_pipes = lost_pipes + pure_pipes
    print(f"Step 5: Calculate the total number of pipes in the cathedral.")
    print(f"Total Pipes = Lost Pipes + Pure Pipes = {lost_pipes} + {pure_pipes} = {total_pipes}.\n")
    
    # Step 6: Answer the question from the poem.
    print("Step 6: Answer the final question.")
    print("How many pipes must the tuner find if half of the lost ones realign?")
    tuner_finds = lost_pipes // 2
    
    # The final print statement shows the equation as requested.
    print(f"The number of pipes for the tuner to find is half of the lost pipes:")
    print(f"{lost_pipes} / 2 = {tuner_finds}")

solve_cathedral_riddle()
<<<140>>>