import math

def solve_seating_arrangement():
    """
    Calculates the number of ways to arrange the group at the circular table.
    """
    
    # Helper for factorial calculation
    def factorial(n):
        return math.factorial(n)

    # --- Part 1: Define group numbers ---
    # Scientists: 12 total (6 male, 6 female). 2 male rowers. 10 others (4 male, 6 female)
    # Mathematicians: 4 total (2 male, 2 female). 1 male rower. 3 others (1 male, 2 female)
    # Ethicists: 2 total (1 male, 1 female)
    # Classicists: 5 total (Cassie is female, 4 others are 2 male, 2 female)

    print("The plan is to create a 'super block' of 16 scientists and mathematicians, and then arrange the 7 other people around it.")
    print("The total number of arrangements is the sum of three cases based on who sits next to this super block.\n")

    # --- Part 2: Calculate possibilities for each case ---

    # Case 1: Both seats next to the super block are filled by ethicists.
    # In this case, the gender of the person at the ends of the super block doesn't matter.
    # Internal SM ways = (Ways to arrange 10 non-rower scientists) * (Ways to arrange 2 rower scientists) * (Ways to arrange 3 non-rower mathematicians)
    num_internal_sm_total = factorial(10) * factorial(2) * factorial(3)
    # External ways = (Ways to arrange 2 ethicists in adjacent seats) * (Ways to arrange 5 classicists in the remaining 5 seats)
    ways_case1 = num_internal_sm_total * factorial(2) * factorial(5)

    # Case 2: Cassie sits next to the Scientist-end, an Ethicist sits next to the Mathematician-end.
    # The scientist at the end of the block must be female (6 women available from the 10 non-rowers).
    # Internal SM ways (S-end female) = (Choose 1 of 6 females for the end) * (Arrange other 9 scientists) * (Arrange 2 rowers) * (Arrange 3 mathematicians)
    num_internal_S_end_F = (6 * factorial(9) * factorial(2)) * factorial(3)
    # External ways = (Place Cassie) * (Place 1 of 2 ethicists) * (Arrange remaining 5 people)
    ways_case2 = num_internal_S_end_F * 2 * factorial(5)

    # Case 3: Cassie sits next to the Mathematician-end, an Ethicist sits next to the Scientist-end.
    # The mathematician at the end must be female (2 women available from the 3 non-rowers).
    # Internal SM ways (M-end female) = (Arrange 10 scientists) * (Arrange 2 rowers) * (Choose 1 of 2 females for end) * (Arrange other 2 mathematicians)
    num_internal_M_end_F = (factorial(10) * factorial(2)) * (2 * factorial(2))
    # External ways = (Place Cassie) * (Place 1 of 2 ethicists) * (Arrange remaining 5 people)
    ways_case3 = num_internal_M_end_F * 2 * factorial(5)
    
    # --- Part 3: Sum the cases and print the results ---
    
    total_ways = ways_case1 + ways_case2 + ways_case3

    print("#--- Calculation Breakdown ---#\n")
    print(f"Subtotal for Case 1 (Ethicist buffers): {num_internal_sm_total} * {factorial(2) * factorial(5)} = {ways_case1}")
    print(f"Subtotal for Case 2 (Cassie by Scientist): {num_internal_S_end_F} * {2 * factorial(5)} = {ways_case2}")
    print(f"Subtotal for Case 3 (Cassie by Mathematician): {num_internal_M_end_F} * {2 * factorial(5)} = {ways_case3}\n")
    
    print("#--- Final Equation ---#\n")
    print(f"Total Ways = (Arrangements for Case 1) + (Arrangements for Case 2) + (Arrangements for Case 3)")
    print(f"Total Ways = {ways_case1} + {ways_case2} + {ways_case3}")
    print(f"Total Ways = {total_ways}")

solve_seating_arrangement()
<<<23652899328000>>>