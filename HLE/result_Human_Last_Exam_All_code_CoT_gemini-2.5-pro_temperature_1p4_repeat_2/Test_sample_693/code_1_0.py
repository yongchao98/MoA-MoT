def calculate_revenue(bids):
    """
    Calculates the total revenue from a single-item divisible auction.

    Args:
        bids: A list of numbers representing the bids.

    Returns:
        The total revenue from the auction.
    """
    if not bids:
        return 0

    # Find the highest bid to determine the winning threshold
    highest_bid = max(bids)
    threshold = 0.5 * highest_bid

    # Identify the winning bids
    winning_bids = [b for b in bids if b >= threshold]

    total_revenue = 0
    # Calculate the payment for each winner
    for winner_bid in winning_bids:
        # Create a list of all other bids
        # We need to handle duplicate bids carefully
        temp_bids = list(bids)
        temp_bids.remove(winner_bid)
        
        if not temp_bids:
            # This case happens if there's only one bidder
            max_other_bid = 0
        else:
            max_other_bid = max(temp_bids)

        # The payment is the critical value needed to win
        payment = 0.5 * max_other_bid
        total_revenue += payment

    return total_revenue

# Bids for the first scenario to calculate x
bids_x = [100, 20, 5]
x = calculate_revenue(bids_x)

# Bids for the second scenario to calculate y
bids_y = [100, 60, 10]
y = calculate_revenue(bids_y)

print(f"For bids {bids_x}, the revenue x is: {int(x)}")
print(f"For bids {bids_y}, the revenue y is: {int(y)}")
print(f"The final answer is (x, y) = ({int(x)}, {int(y)})")
