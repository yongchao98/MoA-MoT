def find_worst_openable_jack():
    """
    This function determines the worst suited Jack to open-raise from the button
    in a standard 100BB NLHE cash game scenario.
    """
    
    # 1. List all suited Jack hands, ordered from best (JTs) to worst (J2s).
    suited_jacks_ranked = ["JTs", "J9s", "J8s", "J7s", "J6s", "J5s", "J4s", "J3s", "J2s"]
    
    # 2. Define a standard GTO button opening range. In a 100BB game, especially
    #    rake-free, this range is very wide and includes all suited Jacks.
    #    We represent the suited Jack portion of the range here.
    button_opening_range = {"JTs", "J9s", "J8s", "J7s", "J6s", "J5s", "J4s", "J3s", "J2s"}
    
    # 3. Find the worst hand that should be opened.
    #    We iterate through our ranked list. The last hand we encounter that is also
    #    in the opening range is by definition the worst one to open.
    worst_jack_to_open = None
    for hand in suited_jacks_ranked:
        if hand in button_opening_range:
            # This hand is a standard open. Since we are iterating from best to worst,
            # we continuously update our answer. The last one will be the correct one.
            worst_jack_to_open = hand
            
    # 4. Print the result.
    if worst_jack_to_open:
        print(f"Based on standard GTO ranges, all suited Jacks are opened from the button.")
        print(f"The hand with the lowest kicker in this category is the 'worst' profitable open.")
        print(f"The worst suited jack you should open from the button is: {worst_jack_to_open}")
    else:
        print("Could not determine the hand based on the provided ranges.")

find_worst_openable_jack()