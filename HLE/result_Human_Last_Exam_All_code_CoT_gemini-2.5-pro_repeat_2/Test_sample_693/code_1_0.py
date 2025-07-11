import collections

def get_auction_outcome(bids_dict):
    """
    Determines the winners and their shares for a given set of bids.
    Args:
        bids_dict: A dictionary of {bidder_id: bid_value}.
    Returns:
        A dictionary of {bidder_id: share} where share is 1/k or 0.
        k: The number of winners.
    """
    if not bids_dict:
        return {}, 0

    bid_values = list(bids_dict.values())
    highest_bid = max(bid_values) if bid_values else 0
    threshold = highest_bid / 2.0

    winners = {
        bidder_id: bid for bidder_id, bid in bids_dict.items() if bid >= threshold
    }
    
    k = len(winners)
    
    shares = collections.defaultdict(float)
    if k > 0:
        for bidder_id in winners:
            shares[bidder_id] = 1.0 / k
            
    return shares, k

def calculate_welfare_for_subset(auction_bids, bidders_to_value):
    """
    Calculates the total welfare for a specific subset of bidders
    based on the outcome of an auction.
    Args:
        auction_bids: Bids to determine the auction outcome.
        bidders_to_value: Bidders whose welfare we want to sum.
    Returns:
        Total welfare for the specified subset.
    """
    shares, k = get_auction_outcome(auction_bids)
    total_welfare = 0
    for bidder_id, bid in bidders_to_value.items():
        # Add the value (bid * share) for the bidder if they are in the subset
        total_welfare += bid * shares.get(bidder_id, 0)
    return total_welfare

def calculate_revenue(bids):
    """
    Calculates the total VCG revenue for a given list of bids.
    """
    # Using dicts to keep track of bidders, assuming unique bids for simplicity
    # or treating them as unique bidders 1, 2, 3...
    initial_bids_dict = {f"bidder_{i+1}": bid for i, bid in enumerate(bids)}
    
    total_revenue = 0
    
    print(f"Calculating revenue for bids: {bids}")
    
    # Calculate payment for each bidder
    for bidder_id, bid_value in initial_bids_dict.items():
        
        # Bids of everyone else
        other_bidders_dict = initial_bids_dict.copy()
        del other_bidders_dict[bidder_id]

        # Case 1: Welfare of others when this bidder participates
        welfare_others_with = calculate_welfare_for_subset(
            initial_bids_dict, other_bidders_dict
        )
        
        # Case 2: Welfare of others when this bidder does NOT participate
        welfare_others_without = calculate_welfare_for_subset(
            other_bidders_dict, other_bidders_dict
        )
        
        # Payment is the "harm" caused to others
        payment = welfare_others_without - welfare_others_with
        
        if payment > 0:
            print(f"  - Payment for bidder with bid {bid_value}: {welfare_others_without:.2f} - {welfare_others_with:.2f} = {payment:.2f}")

        total_revenue += payment
        
    print(f"Total revenue is: {total_revenue:.2f}\n")
    return total_revenue

if __name__ == "__main__":
    # Scenario x
    bids_x = [100, 20, 5]
    x = calculate_revenue(bids_x)

    # Scenario y
    bids_y = [100, 60, 10]
    y = calculate_revenue(bids_y)

    print(f"The final result (x, y) is: ({int(x)}, {int(y)})")
    # Final answer format
    # print(f"<<<({int(x)}, {int(y)})>>>")