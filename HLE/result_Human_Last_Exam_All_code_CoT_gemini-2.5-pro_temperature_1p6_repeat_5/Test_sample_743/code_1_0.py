import itertools

def find_tennis_matchups():
    """
    Finds the largest list of tennis matchups (groups of 4) from 11 players
    such that no two matchups have more than two players in common.
    
    This is solved by finding all 4-player combinations from a set of 11 players
    (labeled 0-10) where the sum of their labels is a multiple of 11.
    This construction guarantees the intersection property.
    """
    players = list(range(11))
    matchups = []

    # Generate all combinations of 4 players from the 11
    all_possible_groups = itertools.combinations(players, 4)

    for group in all_possible_groups:
        # Check if the sum of the player numbers is a multiple of 11
        if sum(group) % 11 == 0:
            matchups.append(group)

    print("Found the following matchups (groups of 4 players):")
    for group in matchups:
        p1, p2, p3, p4 = group
        # Output each number in the final equation
        print(f"  - Group {{{p1}, {p2}, {p3}, {p4}}}: {p1} + {p2} + {p3} + {p4} = {sum(group)}")

    print(f"\nThe largest list of matchups that can be created this way has a size of: {len(matchups)}")

find_tennis_matchups()
<<<22>>>