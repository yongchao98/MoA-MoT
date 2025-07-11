def calculate_social_welfare(bids):
    """
    Calculates allocation and social welfare for a given set of bids.
    A bid represents the bidder's value.
    Returns: a dictionary of {bid: allocation_share} and total_welfare.
    """
    if not bids:
        return {}, 0

    highest_bid = max(bids)
    threshold = highest_bid / 2.0
    
    winners = [b for b in bids if b >= threshold]
    num_winners = len(winners)

    if num_winners == 0:
        return {b: 0 for b in bids}, 0

    share = 1.0 / num_winners
    allocation = {}
    
    # Use a copy to handle duplicate bids correctly
    bids_copy = list(bids)
    for b in bids_copy:
        if b in winners:
            allocation[b] = allocation.get(b, 0) + share
            winners.remove(b)
        else:
            allocation[b] = allocation.get(b, 0) # ensure it has a key

    total_welfare = sum(bid * allocation.get(bid, 0) for bid in bids)
    
    return allocation, total_welfare


def calculate_vcg_revenue(bids):
    """
    Calculates the total VCG revenue for a single-item divisible auction.
    """
    indexed_bids = list(enumerate(bids)) # (original_index, bid_value)
    
    # 1. Calculate allocation and winners with all bidders present
    if not bids:
        print("Revenue = 0")
        return 0

    highest_bid = max(bids)
    threshold = highest_bid / 2.0
    
    winning_bidders = [ib for ib in indexed_bids if ib[1] >= threshold]

    if not winning_bidders:
        print("Revenue = 0")
        return 0

    num_winners = len(winning_bidders)
    share = 1.0 / num_winners

    # 2. Calculate payment for each winner
    total_revenue = 0
    payments = []
    for winner_idx, winner_bid in winning_bidders:
        # a. Calculate welfare of others WITH this winner present
        welfare_others_with_winner = 0
        for idx, bid in indexed_bids:
            if idx != winner_idx:
                # Check if this 'other' bidder is also a winner
                if bid >= threshold:
                    welfare_others_with_winner += bid * share
        
        # b. Calculate welfare of others WITHOUT this winner
        bids_without_winner = [ib[1] for ib in indexed_bids if ib[0] != winner_idx]
        _, welfare_without_winner = calculate_social_welfare(bids_without_winner)

        # c. VCG Payment
        payment = welfare_without_winner - welfare_others_with_winner
        payments.append(payment)
        total_revenue += payment
    
    # Print the detailed calculation
    payment_str = " + ".join([f"{p:.0f}" for p in payments])
    if len(payments) > 1:
        print(f"Revenue = {payment_str} = {total_revenue:.0f}")
    else:
        print(f"Revenue = {total_revenue:.0f}")
        
    return total_revenue

def solve():
    """
    Solves the user's auction problem.
    """
    bids_x = [100, 20, 5]
    bids_y = [100, 60, 10]

    print("Calculating x for bids {}:".format(tuple(bids_x)))
    x = calculate_vcg_revenue(bids_x)

    print("\nCalculating y for bids {}:".format(tuple(bids_y)))
    y = calculate_vcg_revenue(bids_y)

    print(f"\nThe final result for (x, y) is: ({x:.0f}, {y:.0f})")

solve()
<<<20, 80>>>