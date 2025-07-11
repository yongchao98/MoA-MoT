import collections

def _run_auction(bids_dict):
    """
    Runs the auction for a given set of bids.
    
    Args:
        bids_dict: A dictionary of {bidder_index: bid_value}.
        
    Returns:
        A tuple containing:
        - allocation: A dictionary {winner_index: share}.
        - welfare_by_bidder: A dictionary {winner_index: welfare}.
    """
    if not bids_dict:
        return {}, {}

    highest_bid = 0
    for bid in bids_dict.values():
        if bid > highest_bid:
            highest_bid = bid

    threshold = 0.5 * highest_bid
    
    winners = {idx: bid for idx, bid in bids_dict.items() if bid >= threshold}
    
    num_winners = len(winners)
    
    if num_winners == 0:
        return {}, {}
        
    share = 1.0 / num_winners
    
    allocation = {idx: share for idx in winners.keys()}
    welfare_by_bidder = {idx: bid * share for idx, bid in winners.items()}
    
    return allocation, welfare_by_bidder

def calculate_revenue(bids):
    """
    Calculates the total revenue of the truthful auction using VCG payments.
    
    Args:
        bids: A list of bid values.
        
    Returns:
        The total revenue of the auction.
    """
    initial_bids_dict = {i: bids[i] for i in range(len(bids))}
    
    main_allocation, main_welfare_by_bidder = _run_auction(initial_bids_dict)
    
    main_winners = main_allocation.keys()
    
    if not main_winners:
        return 0

    total_revenue = 0
    
    # Calculate payment for each winner
    for winner_idx in main_winners:
        # 1. Welfare of other bidders WITH this winner present
        others_welfare_with_winner = 0
        for other_idx, welfare in main_welfare_by_bidder.items():
            if other_idx != winner_idx:
                others_welfare_with_winner += welfare

        # 2. Welfare of other bidders WITHOUT this winner
        bids_without_winner = initial_bids_dict.copy()
        del bids_without_winner[winner_idx]
        
        _, welfare_dict_without_winner = _run_auction(bids_without_winner)
        total_welfare_without_winner = sum(welfare_dict_without_winner.values())
        
        # 3. Payment is the difference (harm caused)
        payment = total_welfare_without_winner - others_welfare_with_winner
        total_revenue += payment
        
    return total_revenue

# Bids for the first scenario (x)
bids_x = [100, 20, 5]
x = calculate_revenue(bids_x)

# Bids for the second scenario (y)
bids_y = [100, 60, 10]
y = calculate_revenue(bids_y)

# Output the final result in the specified format
# Each number is calculated and printed.
x_val = int(x)
y_val = int(y)
print(f"(x, y) = ({x_val}, {y_val})")