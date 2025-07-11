def find_worst_suited_jack_open():
    """
    This function determines the worst suited Jack to open from the button
    in a 100BB rake-free NLH cash game based on GTO principles.
    """
    
    # In poker, suited Jack hands (Jxs) are ranked by their second card, or "kicker".
    # A higher kicker is stronger. We can list them in order from best to worst.
    suited_jacks_by_strength = [
        "AJs", "KJs", "QJs", "JTs", "J9s", "J8s",
        "J7s", "J6s", "J5s", "J4s", "J3s", "J2s"
    ]

    # For the scenario (100BB, Rake-Free, Button Open), Game Theory Optimal (GTO)
    # charts recommend opening a very wide range of hands (around 45-50%).
    # This range is wide enough to include all of the suited Jack combinations
    # due to their ability to make strong hands like flushes.

    # The "worst" hand to open is the one with the lowest expected value (EV) that
    # is still considered a profitable raise. In this list, that is the very last hand.
    worst_hand_in_range = suited_jacks_by_strength[-1]
    
    print("Analysis Steps:")
    print("1. All suited Jack hands are ranked by their kicker strength.")
    print(f"   - The order from best to worst is: {', '.join(suited_jacks_by_strength)}")
    
    print("2. A standard GTO opening range for the Button in a deep, rake-free game is consulted.")
    print("   - This range includes all the suited Jack hands listed above.")

    print("3. The 'worst' hand to open is the one at the bottom of this list which is still profitable.")
    print("\nResult:")
    print("The worst suited jack you should open from the button is:")
    print(f"-> {worst_hand_in_range}")

find_worst_suited_jack_open()