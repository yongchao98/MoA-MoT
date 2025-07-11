def get_auction_outcome(bids):
    """
    Determines the winners and their allocation based on the auction rule.
    Args:
        bids (dict): A dictionary of {bidder_id: bid_value}.
    Returns:
        dict: A dictionary of {winner_id: allocation_share}.
    """
    if not bids:
        return {}

    max_bid = max(bids.values())
    threshold = max_bid / 2.0

    winners = {bidder_id: bid for bidder_id, bid in bids.items() if bid >= threshold}

    if not winners:
        return {}

    num_winners = len(winners)
    share = 1.0 / num_winners

    allocation = {winner_id: share for winner_id in winners}
    return allocation

def calculate_vcg_revenue(initial_bids):
    """
    Calculates the total auction revenue using the VCG payment rule.
    Args:
        initial_bids (list): A list of bid values.
    Returns:
        float: The total revenue of the auction.
    """
    # Use dictionary with original indices to track bidders
    bids_dict = {i: bid for i, bid in enumerate(initial_bids)}

    # 1. Determine the outcome with all bidders
    initial_allocation = get_auction_outcome(bids_dict)
    winners = initial_allocation.keys()
    
    total_revenue = 0

    # 2. For each winner, calculate their VCG payment
    for winner_id in winners:
        # Calculate welfare of others WITH this winner
        welfare_with_winner = 0
        for bidder_id, bid in bids_dict.items():
            if bidder_id != winner_id:
                welfare_with_winner += bid * initial_allocation.get(bidder_id, 0)

        # Calculate welfare of others WITHOUT this winner
        bids_without_winner = {bidder_id: bid for bidder_id, bid in bids_dict.items() if bidder_id != winner_id}
        allocation_without_winner = get_auction_outcome(bids_without_winner)
        
        welfare_without_winner = 0
        for bidder_id, bid in bids_without_winner.items():
            welfare_without_winner += bid * allocation_without_winner.get(bidder_id, 0)
            
        # Payment is the "harm" caused to others
        payment = welfare_without_winner - welfare_with_winner
        total_revenue += payment
        
    return total_revenue

# Bids for the two scenarios
bids_x = [100, 20, 5]
bids_y = [100, 60, 10]

# Calculate revenues
x = calculate_vcg_revenue(bids_x)
y = calculate_vcg_revenue(bids_y)

# Print the final result in the specified format
print(f"(x, y) = ({int(x)}, {int(y)})")
<<< (20, 80) >>>