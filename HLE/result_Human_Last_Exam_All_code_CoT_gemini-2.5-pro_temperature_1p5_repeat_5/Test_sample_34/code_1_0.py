def find_worst_suited_jack_to_open():
    """
    This script determines the worst suited Jack to open-raise from the button
    in a 100BB cash game based on standard Game Theory Optimal (GTO) strategy.
    """

    # In GTO poker theory, the button's opening range is wide and includes
    # all suited Jack combinations because they are all profitable raises.
    gto_opening_range_suited_jacks = [
        "JAs", "JKs", "JQs", "JTs", "J9s", "J8s",
        "J7s", "J6s", "J5s", "J4s", "J3s", "J2s"
    ]

    # To find the "worst" hand, we find the one with the lowest-ranked kicker.
    # We assign numerical values to each card to perform this comparison.
    card_ranks = {
        'A': 14, 'K': 13, 'Q': 12, 'T': 10, '9': 9,
        '8': 8, '7': 7, '6': 6, '5': 5, '4': 4, '3': 3, '2': 2
    }

    worst_hand_found = ""
    # We start with a very high number to ensure the first hand's rank is lower.
    lowest_kicker_rank = float('inf')

    print("Analyzing the kicker strength of all suited Jacks in the button's opening range:")
    print("--------------------------------------------------------------------------")

    # Iterate through the list of hands to find the one with the weakest kicker.
    for hand in gto_opening_range_suited_jacks:
        # The kicker is the second character in the hand's name (e.g., 'A' in "JAs").
        kicker_character = hand[1]
        kicker_rank_value = card_ranks[kicker_character]

        # Output the hand and its kicker's numerical rank for comparison.
        print(f"Hand: {hand}, Kicker: '{kicker_character}', Kicker Rank: {kicker_rank_value}")

        # Check if this kicker's rank is the lowest we have seen so far.
        if kicker_rank_value < lowest_kicker_rank:
            lowest_kicker_rank = kicker_rank_value
            worst_hand_found = hand

    print("--------------------------------------------------------------------------")
    print(f"The hand with the minimum kicker rank ({lowest_kicker_rank}) is the 'worst' suited Jack to open.")
    print(f"\nThe worst suited jack you should open from the button is: {worst_hand_found}")


if __name__ == '__main__':
    find_worst_suited_jack_to_open()