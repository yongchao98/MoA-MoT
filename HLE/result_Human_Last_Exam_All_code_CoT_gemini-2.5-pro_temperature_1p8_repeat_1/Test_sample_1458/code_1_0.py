def solve_cathedral_echo():
    """
    This function solves the riddle from "The Cathedral's Echo".

    Let L be the number of lost pipes (out of tune).
    Let T be the total number of pipes.
    The number of pure (in tune) pipes is 200.

    From the poem, we can form two equations:
    1) L = T / 3  ("One-third of pipes fell out of tune")
    2) T = L + 200 ("two hundred pipes still sing pure")

    Substitute (2) into (1):
    L = (L + 200) / 3
    3 * L = L + 200
    2 * L = 200
    L = 100

    The other fractions (3/7, 1/4) are considered poetic red herrings,
    as they are not mathematically consistent with the primary facts (i.e., 100 is not divisible by 7).
    """

    # We solve the equation 2 * L = 200 to find the number of lost pipes.
    lost_pipes = 200 // 2

    # The question asks for the number of pipes the tuner must find, which is half of the lost pipes.
    tuner_finds = lost_pipes // 2
    
    print(f"First, we determine the number of lost pipes. The relationship is T = L + 200 and L = T / 3. This resolves to 2 * L = 200, so the number of lost pipes is {lost_pipes}.")
    print("The final question asks for half of the lost pipes that the tuner must find.")
    print(f"The final calculation is:")
    print(f"{lost_pipes} / 2 = {tuner_finds}")

solve_cathedral_echo()