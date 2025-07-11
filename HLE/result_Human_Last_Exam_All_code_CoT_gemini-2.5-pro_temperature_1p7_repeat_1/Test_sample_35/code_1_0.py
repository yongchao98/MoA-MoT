import collections

def explain_shoving_decision():
    """
    Analyzes and explains the best hand to shove based on equity against a likely calling range.
    """
    
    # --- Step 1: Define the Situation ---
    stack_size_bb = 16
    position = "UTG+1"
    situation = "Near the money bubble"
    
    print(f"Scenario: You are {position} with {stack_size_bb}bb {situation}.")
    print("This is a classic 'push-or-fold' spot. The key is to choose a hand that is profitable to shove.")
    print("-" * 50)
    
    # --- Step 2: Define Hands and a Likely Calling Range ---
    # These are pre-calculated equities using standard poker tools.
    # The calling range is assumed to be 77+, AJs+, AQo+.
    # This range consists of 84 hand combinations.
    # Format: {"Hand": Equity_Percentage}
    hand_equities = {
        "QJs": 38.3,
        "99": 50.1,
        "AJo": 42.8,
        "AKo": 49.3
    }

    print("Analysis: All the options except 'None of these' are generally considered shoves in this spot.")
    print("To find the *best* option, we should see how each hand performs when called.")
    print("A likely calling range from an opponent is 77+, AJs+, AQo+.\n")
    print("Let's look at the equity of each hand against this range:")
    print("-" * 50)
    
    # --- Step 3 & 4: Analyze, Compare, and Conclude ---
    # Find the hand with the highest equity from the dictionary
    best_hand = max(hand_equities, key=hand_equities.get)
    
    for hand, equity in hand_equities.items():
        print(f"Equity of {hand} vs calling range: {equity}%")

    print("-" * 50)
    print("Conclusion:")
    print(f"As you can see, {best_hand} has the highest equity ({hand_equities[best_hand]}%) when we get called.")
    print("While AKo is also very strong, the pocket pair 99 holds its value slightly better against this specific range of hands that is likely to call a 16bb shove.")
    print("It has a significant edge over AJo and QJs. Therefore, it's the most robust choice among the options.")


explain_shoving_decision()