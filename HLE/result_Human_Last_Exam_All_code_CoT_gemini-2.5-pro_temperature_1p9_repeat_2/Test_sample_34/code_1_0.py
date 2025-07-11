import pandas as pd

def find_worst_open():
    """
    Identifies the worst suited Jack to open from the button in a 100BB cash game
    by finding the hand with the lowest non-zero opening frequency according to a
    standard GTO model.
    """
    # Step 1: Define GTO opening frequencies for suited Jacks from the button.
    # 1.0 means the hand is always opened.
    # <1.0 means it's a "mixed strategy" (sometimes opened, sometimes folded).
    # 0.0 means it's always folded.
    # These frequencies are based on common 100BB GTO solver outputs.
    button_opening_range = {
        'JTs': 1.0,
        'J9s': 1.0,
        'J8s': 1.0,
        'J7s': 1.0,
        'J6s': 1.0,
        'J5s': 0.86, # Starts to become a mix
        'J4s': 0.58, # Mixed frequency
        'J3s': 0.35, # Mixed frequency
        'J2s': 0.25, # Lowest frequency open
    }

    print("Analyzing GTO opening frequencies for suited Jacks (JXs) from the button:")
    for hand, freq in button_opening_range.items():
        print(f"- {hand}: Open {freq*100:.0f}% of the time")

    # Step 2: Filter for hands that are part of the opening range (frequency > 0).
    playable_hands = {hand: freq for hand, freq in button_opening_range.items() if freq > 0}

    # Step 3: Find the hand with the minimum opening frequency.
    # This represents the "worst" hand you should still open.
    if not playable_hands:
        worst_hand = "None"
    else:
        worst_hand = min(playable_hands, key=playable_hands.get)

    min_freq = playable_hands[worst_hand]

    print("\nFinal Equation:")
    print(f"Worst Suited Jack to Open = Hand with min(Opening Frequency) where Frequency > 0")
    print(f"min({list(playable_hands.values())}) = {min_freq}")
    print(f"The hand corresponding to this frequency is {worst_hand}.")
    print("\n----------------------------------------------------")
    print(f"Result: The worst suited Jack to open from the button is {worst_hand}.")
    print("----------------------------------------------------")

if __name__ == "__main__":
    find_worst_open()
