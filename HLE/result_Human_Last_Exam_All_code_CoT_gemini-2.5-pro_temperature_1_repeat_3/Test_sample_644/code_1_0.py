# This script explains and simulates the ritual to reverse the bad luck
# from climbing the Leaning Tower of Pisa before graduation.

def solve_pisa_superstition():
    """
    Explains and simulates the superstitious fix for students in Pisa.
    """
    print("You have made the mistake of climbing the Leaning Tower of Pisa before graduating.")
    print("According to student superstition, this means you will not graduate!")
    print("\nHowever, there is a traditional way to fix this and reverse the bad luck.")
    print("The solution involves the large, round Baptistery located next to the Tower.")
    print("\nTHE FIX: You must go around the Baptistery 100 times.")
    print("\nLet's simulate completing this ritual to break the curse.")
    print("-" * 30)

    total_laps = 100
    # The 'equation' to solve is completing 100 laps.
    # We will now output each number in this sequence.
    for lap_number in range(1, total_laps + 1):
        print(f"Lap {lap_number} of {total_laps} complete...")

    print("-" * 30)
    print("\nCongratulations! You have completed all 100 laps around the Baptistery.")
    print("According to the legend, the curse is now broken. Good luck with your graduation!")

# Run the simulation to find the answer.
solve_pisa_superstition()