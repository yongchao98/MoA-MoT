def calculate_auction_revenue(bids):
    """
    Calculates the revenue of the described divisible item auction.

    The allocation rule is that the item is divided among all bidders whose
    bid is at least 1/2 of the highest bid. For the auction to be truthful,
    the revenue is the highest bid among those who did not win a share.

    Args:
        bids: A list of numbers representing the bids.

    Returns:
        The total revenue of the auction.
    """
    if not bids:
        return 0

    highest_bid = max(bids)
    winning_threshold = highest_bid / 2.0

    winning_bids = []
    losing_bids = []

    for bid in bids:
        if bid >= winning_threshold:
            winning_bids.append(bid)
        else:
            losing_bids.append(bid)

    if not losing_bids:
        # If everyone wins, the clearing price is 0.
        revenue = 0
    else:
        # Revenue is the highest losing bid.
        revenue = max(losing_bids)

    return revenue

# Scenario 1 for x
bids_x = [100, 20, 5]
x = calculate_auction_revenue(bids_x)

# Scenario 2 for y
bids_y = [100, 60, 10]
y = calculate_auction_revenue(bids_y)

# Print the results step-by-step
print(f"For bids {tuple(bids_x)}, the revenue x is: {int(x)}")
print(f"For bids {tuple(bids_y)}, the revenue y is: {int(y)}")
print(f"The result (x, y) is: ({int(x)}, {int(y)})")
