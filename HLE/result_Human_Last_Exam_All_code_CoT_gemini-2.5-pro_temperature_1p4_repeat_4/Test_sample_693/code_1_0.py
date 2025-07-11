import collections

def run_auction(bids_dict):
    """
    Runs the auction based on the described allocation rule.
    
    Args:
        bids_dict: A dictionary of {bidder_id: bid_value}.
        
    Returns:
        A tuple containing:
        - winners: A list of bidder_ids who won.
        - allocation: A dictionary of {bidder_id: share} for all bidders.
    """
    if not bids_dict:
        return [], collections.defaultdict(float)
        
    highest_bid = max(bids_dict.values())
    threshold = highest_bid / 2.0
    
    winners = [bidder_id for bidder_id, bid in bids_dict.items() if bid >= threshold]
    
    allocation = collections.defaultdict(float)
    if not winners:
        return [], allocation
        
    num_winners = len(winners)
    share = 1.0 / num_winners
    
    for winner_id in winners:
        allocation[winner_id] = share
        
    return winners, allocation

def calculate_vcg_revenue(bids):
    """
    Calculates the total revenue of the auction using VCG payments.
    
    Args:
        bids: A list of bid values.
        
    Returns:
        The total revenue.
    """
    # Use bids as valuations, assuming truthfulness
    valuations = {i: bid for i, bid in enumerate(bids)}
    
    # Run the auction with all bidders
    winners, allocation_with_all = run_auction(valuations)
    
    if not winners:
        return 0
        
    total_revenue = 0
    
    print(f"Bids: {bids}")
    print(f"Highest bid: {max(bids)}. Threshold: {max(bids)/2.0}")
    print(f"Winners are bidders with bids: {[bids[w] for w in winners]}")
    if winners:
      print(f"Each of the {len(winners)} winners gets 1/{len(winners)} of the item.")
    print("-" * 20)
    
    # Calculate payment for each winner
    for winner_id in winners:
        # 1. Calculate welfare of others WITH the winner present
        welfare_of_others_with_winner = 0
        for bidder_id, share in allocation_with_all.items():
            if bidder_id != winner_id:
                welfare_of_others_with_winner += valuations[bidder_id] * share
        
        # 2. Calculate welfare of others WITHOUT the winner present
        bids_without_winner = valuations.copy()
        del bids_without_winner[winner_id]
        
        _, allocation_without_winner = run_auction(bids_without_winner)
        
        welfare_of_others_without_winner = 0
        for bidder_id, share in allocation_without_winner.items():
            # Note: bidder_id here are from the dict without the current winner
            welfare_of_others_without_winner += valuations[bidder_id] * share
            
        # 3. Payment is the difference (harm caused)
        payment = welfare_of_others_without_winner - welfare_of_others_with_winner
        total_revenue += payment
        
        print(f"Payment for winner with bid {valuations[winner_id]}:")
        print(f"  Welfare of others without this winner = {welfare_of_others_without_winner:.2f}")
        print(f"  Welfare of others with this winner = {welfare_of_others_with_winner:.2f}")
        print(f"  Payment = {welfare_of_others_without_winner:.2f} - {welfare_of_others_with_winner:.2f} = {payment:.2f}")
        print("-" * 20)

    return total_revenue

def solve_problem():
    """
    Solves the problem by calculating revenues x and y.
    """
    bids_x = [100, 20, 5]
    bids_y = [100, 60, 10]
    
    print("--- Calculating revenue x ---")
    x = calculate_vcg_revenue(bids_x)
    print(f"\nTotal Revenue x = {x:.2f}")
    
    print("\n" + "="*30 + "\n")
    
    print("--- Calculating revenue y ---")
    y = calculate_vcg_revenue(bids_y)
    print(f"\nTotal Revenue y = {y:.2f}")
    
    print("\n" + "="*30 + "\n")
    
    print(f"The final result for (x, y) is ({x:.0f}, {y:.0f})")

solve_problem()
<<<20, 80>>>