def calculate_shove_ev(hand_name, equity_vs_calling_range):
    """
    Calculates and prints the Expected Value (EV) of shoving a specific hand.
    
    The EV formula is:
    EV = (P_fold * pot_before_shove) + (P_call * (our_equity * total_pot_when_called - amount_risked))
    """
    
    # --- 1. Define Scenario Parameters ---
    our_stack = 16.0  # Our stack in big blinds
    pot_before_shove = 2.5  # SB (0.5) + BB (1.0) + Antes (1.0, estimate)
    p_fold = 0.70  # Probability all players fold (high due to bubble pressure)
    p_call = 1.0 - p_fold # Probability we get called
    
    # --- 2. Calculate Values for the EV Equation ---
    amount_risked = our_stack
    total_pot_when_called = pot_before_shove + our_stack + our_stack # Pot + our shove + caller's stack
    
    # --- 3. Calculate the two components of EV ---
    fold_equity_value = p_fold * pot_before_shove
    call_equity_value = p_call * (equity_vs_calling_range * total_pot_when_called - amount_risked)
    
    # --- 4. Calculate Final EV ---
    total_ev = fold_equity_value + call_equity_value
    
    # --- 5. Print the breakdown ---
    print(f"--- Analysis for {hand_name} ---")
    print(f"Equity vs. a tight calling range (99+, AQs+, AKo): {equity_vs_calling_range:.2%}")
    # The prompt requires printing each number in the final equation.
    print("EV Equation: (P_fold * Pot) + (P_call * (Equity * Total_Pot - Risk))")
    print(f"EV = ({p_fold:.2f} * {pot_before_shove}) + ({p_call:.2f} * ({equity_vs_calling_range:.2f} * {total_pot_when_called} - {amount_risked}))")
    print(f"EV = {fold_equity_value:.3f} + ({p_call:.2f} * ({equity_vs_calling_range * total_pot_when_called:.2f} - {amount_risked}))")
    print(f"EV = {fold_equity_value:.3f} + ({call_equity_value:.3f})")
    print(f"Final EV for {hand_name}: {total_ev:.3f} big blinds\n")
    return total_ev

def main():
    """
    Analyzes the poker hand choices to find the most profitable shove.
    """
    print("Poker Hand Analysis: 16bb shove from UTG+1 near the money bubble.\n")

    # Equities are calculated against a tight calling range of {99+, AQs+, AKo}
    # These are standard equity values from poker software.
    hand_equities = {
        "A. QJs": 0.30,
        "C. 99": 0.36,  # Equity is lower because when called by 99+, it's often against a better pair
        "D. AJo": 0.26,
        "E. AKo": 0.44
    }

    results = {}
    for hand, equity in hand_equities.items():
        ev = calculate_shove_ev(hand, equity)
        results[hand] = ev
        
    # --- Conclusion ---
    best_hand = max(results, key=results.get)
    print("--- Conclusion ---")
    print("Based on the Expected Value (EV) calculations:")
    print(" - QJs is barely profitable but very risky.")
    print(" - AJo is unprofitable (negative EV) and should be folded.")
    print(" - 99 is a clearly profitable shove.")
    print(" - AKo is also clearly profitable and yields the highest return on average.")
    print(f"\nThe best hand to jam is the one with the highest EV: {best_hand}")

if __name__ == "__main__":
    main()