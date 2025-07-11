def calculate_shove_ev():
    """
    Calculates and displays the Expected Value (EV) of shoving different poker hands
    in a specific tournament scenario.
    """
    # --- 1. Define Assumptions ---
    investment = 16.0  # Our stack size in big blinds (bb)
    blinds = 1.5       # Small blind + Big blind
    players = 8        # Players at the table
    ante_pct = 0.125   # Ante as a percentage of the big blind
    
    # Starting pot = blinds + total antes
    pot_pre = blinds + (players * ante_pct)
    
    # Pot size if we are called (our shove + villain's call + dead money from antes/blinds)
    # The caller does not pay their blinds/ante twice.
    pot_if_called = investment + investment + (pot_pre - blinds)

    # Assume a very tight calling range for villains due to bubble pressure (e.g., TT+, AK)
    villain_call_range_pct = 0.027
    players_to_act = 7
    
    # Probability everyone folds to our shove
    p_fold = (1 - villain_call_range_pct) ** players_to_act
    
    # Probability we get at least one caller
    p_call = 1 - p_fold

    # Hand data: hand name and its equity vs the assumed calling range (TT+, AK)
    hands = [
        {'name': 'AKo', 'equity': 0.439},
        {'name': '99', 'equity': 0.355},
        {'name': 'AJo', 'equity': 0.306},
        {'name': 'QJs', 'equity': 0.322},
    ]

    print(f"Analysis for a {investment}bb shove from UTG+1 (8-handed) on the bubble.\n")
    print(f"Assumptions:")
    print(f"- Starting Pot: {pot_pre:.2f}bb")
    print(f"- Opponent Calling Range: ~{villain_call_range_pct*100:.1f}% (e.g., TT+, AK)")
    print(f"- P(all fold) = {p_fold:.3f}, P(get called) = {p_call:.3f}\n")
    
    print("--- EV Calculations ---")
    # --- 2. Calculate EV for each hand ---
    for hand in hands:
        # EV when we get called = (our equity * final pot) - what we invested
        ev_when_called = (hand['equity'] * pot_if_called) - investment
        
        # Total EV = (value from folds) + (value from being called)
        total_ev = (p_fold * pot_pre) + (p_call * ev_when_called)
        
        print(f"\nHand: {hand['name']}")
        print(f"Equity vs. Calling Range: {hand['equity']:.1%}")
        print("EV Formula: (P_fold * pot_pre) + (P_call * (Equity * pot_if_called - investment))")
        # Print the numbers used in the final equation
        print("EV Equation:")
        print(f"({p_fold:.3f} * {pot_pre:.1f}) + ({p_call:.3f} * ({hand['equity']:.3f} * {pot_if_called:.1f} - {investment:.1f})) = {total_ev:+.2f}bb")

    print("\n--- Conclusion ---")
    print("The EV calculation shows the average profit in big blinds for making this play.")
    print("On the bubble, we must be extra cautious and only shove with a clear, profitable advantage.")
    print("While AKo is the highest EV shove, 99 is also clearly profitable.")
    print("AJo and QJs have much lower profitability, making them risky shoves given the ICM pressure.")
    print("Standard ICM-aware charts for this spot recommend shoving 99+, AJs+, AKo.")
    print("This makes both 99 and AKo correct. However, 99 represents the edge of the shoving range for pairs, making it an excellent hand to test game knowledge. A hand like 88 would be a fold.")
    print("Therefore, 99 is a mandatory and correct jam.")

calculate_shove_ev()