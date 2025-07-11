def calculate_auction_revenue(bids, case_name):
    """
    Calculates the total revenue for a divisible item auction based on the given rules.

    Args:
        bids (list): A list of the bids from all bidders.
        case_name (str): The name of the case ('x' or 'y') for printing.

    Returns:
        float: The total revenue from the auction.
    """
    print(f"--- Calculating revenue {case_name} for bids {bids} ---")
    
    if not bids:
        print("No bids provided. Revenue is 0.")
        return 0.0

    highest_bid = max(bids)
    winning_threshold = 0.5 * highest_bid
    print(f"Highest bid is {highest_bid}, so the winning threshold is {winning_threshold}.")

    # Find the winners
    winners = []
    for i, bid in enumerate(bids):
        if bid >= winning_threshold:
            winners.append((i, bid))

    if not winners:
        print("There are no winners. Total revenue is 0.")
        return 0.0

    print(f"The winners are bidders with bids: {[w[1] for w in winners]}.")
    
    total_revenue = 0
    payment_terms = []
    
    # Calculate payment for each winner
    for i, bid in winners:
        # Bids of all other participants
        other_bids = [b for j, b in enumerate(bids) if i != j]
        
        if not other_bids:
            # If there is only one bidder, they win but pay 0
            payment = 0.0
        else:
            highest_other_bid = max(other_bids)
            # Payment is the critical price needed to win
            payment = 0.5 * highest_other_bid
        
        print(f"Payment for bidder with bid {bid}: 0.5 * max({other_bids}) = {payment}")
        total_revenue += payment
        payment_terms.append(str(payment))

    if len(payment_terms) > 1:
      print(f"Total revenue {case_name} = {' + '.join(payment_terms)} = {total_revenue}")
    else:
      print(f"Total revenue {case_name} = {total_revenue}")
      
    return total_revenue

# Bids for case x
bids_x = [100, 20, 5]
x = calculate_auction_revenue(bids_x, 'x')

print("\n" + "="*30 + "\n")

# Bids for case y
bids_y = [100, 60, 10]
y = calculate_auction_revenue(bids_y, 'y')

print("\n" + "="*30 + "\n")
print(f"The final result for (x, y) is ({x}, {y}).")