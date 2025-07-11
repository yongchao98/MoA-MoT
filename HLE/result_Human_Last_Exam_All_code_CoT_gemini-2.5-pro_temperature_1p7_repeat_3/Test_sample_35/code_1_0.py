import sys

def get_rank(card):
    """Returns the numerical rank of a card."""
    if card.isdigit():
        return int(card)
    return {'T': 10, 'J': 11, 'Q': 12, 'K': 13, 'A': 14}[card]

def is_in_range(hand, hand_range):
    """
    Checks if a hand is within a defined GTO range.
    Range notation examples: '77+' (pairs), 'A9s+' (suited aces), 'KQo' (specific offsuit).
    """
    card1_rank = get_rank(hand[0])
    card2_rank = get_rank(hand[1])
    is_suited = len(hand) == 3 and hand[2] == 's'
    is_pair = card1_rank == card2_rank

    for r_hand in hand_range:
        # Case 1: Pocket Pair Range (e.g., 'TT+')
        if r_hand.endswith('+') and r_hand[0] == r_hand[1]:
            if is_pair:
                min_rank = get_rank(r_hand[0])
                if card1_rank >= min_rank:
                    return True
        # Case 2: Suited High Card Range (e.g., 'AJs+')
        elif r_hand.endswith('s+'):
            if is_suited and card1_rank == get_rank(r_hand[0]):
                min_rank_kicker = get_rank(r_hand[1])
                if card2_rank >= min_rank_kicker:
                    return True
        # Case 3: Offsuit High Card Range (e.g., 'AQo+')
        elif r_hand.endswith('o+'):
            if not is_suited and card1_rank == get_rank(r_hand[0]):
                min_rank_kicker = get_rank(r_hand[1])
                if card2_rank >= min_rank_kicker:
                    return True
        # Case 4: Specific Hand (e.g., AKs, KQo)
        else:
            if hand.lower() == r_hand.lower():
                return True
    return False

def main():
    """
    Main function to analyze the poker hand choices.
    """
    # Define a tight shoving range for 16bb UTG+1 on the bubble.
    # This range accounts for heavy ICM pressure.
    # TT+ means TT, JJ, QQ, KK, AA.
    # AKs is Ace-King suited.
    # AKo is Ace-King offsuit.
    shoving_range = ['TT+', 'AKs', 'AKo']

    # The user's hand options
    options = {
        'A': 'QJs',
        'C': '99',
        'D': 'AJo',
        'E': 'AKo'
    }

    print("Analyzing shove/fold decisions for 16bb UTG+1 on the money bubble.")
    print(f"Using a tight GTO Range: {', '.join(shoving_range)}\n")

    correct_option = None
    for option, hand in options.items():
        if is_in_range(hand, shoving_range):
            action = "SHOVE"
            correct_option = option
        else:
            action = "FOLD"
        
        # We need to print each part of the final equation to be helpful
        # The equation here is the logical check of each option
        print(f"Hand: {hand} ({option}) -> Action: {action}")
        
    # The problem implies a single best answer. 
    # Our tight, bubble-adjusted range yields one correct shove.
    # The following line is not part of the required output format but helps clarify.
    # print(f"\nBased on this analysis, the correct hand to jam is {options.get(correct_option)}.")


if __name__ == "__main__":
    main()
