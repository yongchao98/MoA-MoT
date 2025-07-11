import math

def get_outcome(bids):
    """
    Calculates the outcome of the auction for a given set of bids.
    Returns the list of winning bids and the total social welfare.
    We assume bids are equivalent to valuations for the whole item.
    """
    if not bids:
        return [], 0

    highest_bid = max(bids)
    threshold = 0.5 * highest_bid
    
    winners = sorted([bid for bid in bids if bid >= threshold], reverse=True)
    
    if not winners:
        return [], 0
        
    k = len(winners)
    # Social welfare is the sum of the valuations of the winners for their share.
    # A winner with valuation 'v' for the whole item gets a share 1/k, worth v/k.
    welfare = sum(winner / k for winner in winners)
    
    return winners, welfare

def calculate_vcg_revenue(bids, scenario_name):
    """
    Calculates the total revenue using VCG payments for a given scenario.
    """
    print(f"--- Calculating revenue '{scenario_name}' for bids: {bids} ---")
    
    # 1. Calculate the outcome with all bidders present
    original_winners, original_welfare = get_outcome(bids)
    
    if not original_winners:
        print("No winners in this scenario. Revenue is 0.")
        return 0

    print(f"With all bidders, the highest bid is {max(bids)} and the winning threshold is 0.5 * {max(bids)} = {0.5 * max(bids):.2f}.")
    print(f"The winners are {original_winners}. The total welfare is {original_welfare:.2f}.")
    
    total_revenue = 0
    
    # 2. For each winner, calculate their VCG payment
    for winner in original_winners:
        # Create a copy of the bids to modify
        bids_copy = list(bids)
        
        # Create a list of bids without the current winner
        bids_without_winner = list(bids_copy)
        bids_without_winner.remove(winner)
        
        # Calculate welfare if this winner had not participated
        _, welfare_without_winner = get_outcome(bids_without_winner)
        
        # Calculate the welfare of all *other* agents in the original outcome
        value_winner_gets = winner / len(original_winners)
        welfare_of_others = original_welfare - value_winner_gets
        
        # VCG Payment = (Welfare if winner is absent) - (Welfare of others if winner is present)
        payment = welfare_without_winner - welfare_of_others
        
        print(f"\nCalculating payment for bidder {winner}:")
        print(f"  Welfare if bidder {winner} is absent: {welfare_without_winner:.2f}")
        print(f"  Welfare of others if bidder {winner} is present: {original_welfare:.2f} (total) - {value_winner_gets:.2f} (value to {winner}) = {welfare_of_others:.2f}")
        print(f"  Payment equation: {welfare_without_winner:.2f} - {welfare_of_others:.2f} = {payment:.2f}")
        
        total_revenue += payment
        
    # Using math.ceil to handle potential floating point artifacts and ensure clean integer boundaries
    final_revenue = math.ceil(total_revenue - 1e-9)
    print(f"\nTotal revenue ({scenario_name}) = {final_revenue}")
    print("-" * 40)
    return final_revenue

def solve_auction():
    bids_x = [100, 20, 5]
    bids_y = [100, 60, 10]
    
    x = calculate_vcg_revenue(bids_x, 'x')
    y = calculate_vcg_revenue(bids_y, 'y')
    
    print(f"The final result for (x, y) is: ({x}, {y})")

solve_auction()