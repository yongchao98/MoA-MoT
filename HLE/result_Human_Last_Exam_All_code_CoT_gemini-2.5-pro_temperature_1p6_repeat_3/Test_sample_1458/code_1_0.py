import math

def solve_cathedral_echo():
    """
    Solves the Cathedral's Echo riddle by finding the total number of pipes
    and then calculating the number of pipes the tuner must find.
    """
    # From the poem, the number of pipes that "lost their perfect pitch" (L)
    # must be divisible by both 7 and 4.
    # L = T / 3, so L must be a multiple of lcm(7, 4) = 28.
    # This means T must be a multiple of 3 * 28 = 84.
    base_multiple = 3 * (7 * 4 // math.gcd(7, 4))
    
    total_pipes = 0
    pipes_in_tune_target = 200

    # We test multiples of 84 to find the correct total number of pipes (T).
    # We assume the two groups of changed pipes ('out of tune' and 'caught the moon') overlap.
    # Tuned Pipes = T - (|L| + |M| - |L intersect M|)
    # 200 = T - (T/3 + 2*T/5 - |L intersect M|)
    # 200 = T - 11*T/15 + |L intersect M|
    # |L intersect M| = 200 - (4*T/15)
    
    # We loop to find the first 'T' that provides a valid, whole-number intersection.
    t_candidate = base_multiple
    while True:
        lost_pitch_pipes = t_candidate / 3
        moon_pipes = (2 * t_candidate) / 5

        # The intersection must be a whole number, positive, and smaller than both sets.
        intersection = pipes_in_tune_target - (4 * t_candidate / 15)
        
        if (intersection > 0 and 
            intersection == int(intersection) and 
            intersection <= lost_pitch_pipes and
            intersection <= moon_pipes):
            total_pipes = t_candidate
            break
        
        t_candidate += base_multiple

    # "The lost" refers to the T/3 pipes that "lost their perfect pitch".
    pipes_to_fix = (total_pipes / 3) / 2

    # Outputting the numbers from the equation as requested
    print("Solving the Cathedral's Echo:")
    print(f"The total number of pipes in the organ (T) must be a multiple of 84.")
    print(f"By checking these multiples, we find T = {int(total_pipes)}.")
    print("\nEquation for the final answer:")
    print(f"Number of pipes that lost their pitch = Total Pipes / 3")
    print(f"                                   {int(total_pipes / 3)} = {int(total_pipes)} / 3")
    print(f"\nNumber the tuner must find = Pipes that lost their pitch / 2")
    print(f"Final Answer                 = {int(total_pipes / 3)} / 2 = {int(pipes_to_fix)}")

solve_cathedral_echo()
<<<70>>>