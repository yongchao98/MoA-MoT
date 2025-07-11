def calculate_auction_revenue(bids, case_name):
    """
    Calculates the total revenue of a divisible item auction based on a truthful payment rule.
    
    The allocation rule is that the item is divided equally among all bidders
    whose bid is at least 1/2 of the highest bid.

    The payment rule for a truthful auction makes each winner `i` pay
    the lowest bid they could have made to win, which is 0.5 * max(other bids).
    """
    print(f"--- Calculating Revenue '{case_name}' for Bids: {bids} ---")
    
    if not bids or len(bids) < 2:
        print("Revenue is 0 as there are not enough bidders.")
        print("-" * 30)
        return 0

    highest_bid = max(bids)
    print(f"The highest bid (b_max) is: {highest_bid}")

    winning_threshold = 0.5 * highest_bid
    print(f"The winning threshold (0.5 * b_max) is: {winning_threshold:.2f}")

    # Identify winners by their original index and bid value
    winners = []
    for index, bid in enumerate(bids):
        if bid >= winning_threshold:
            winners.append({'index': index, 'bid': bid})

    if not winners:
        print("There are no winners, so the revenue is 0.")
        print("-" * 30)
        return 0
    
    winning_bids = [w['bid'] for w in winners]
    print(f"Winning Bids are: {winning_bids}")
    
    total_revenue = 0
    payment_equations = []

    print("\nCalculating payments for each winner:")
    for winner in winners:
        # Get bids from all other bidders
        other_bids = [b for i, b in enumerate(bids) if i != winner['index']]
        highest_other_bid = max(other_bids)
        
        # The payment is 0.5 * highest bid among others
        payment = 0.5 * highest_other_bid
        
        print(f"  For winner with bid {winner['bid']}:")
        print(f"    - Highest bid among others is {highest_other_bid}.")
        print(f"    - Payment = 0.5 * {highest_other_bid} = {payment:.2f}")
        
        total_revenue += payment
        payment_equations.append(f"{payment:.2f}")
    
    print("\nCalculating Total Revenue:")
    if len(payment_equations) > 1:
        final_equation = f"Total Revenue ('{case_name}') = {' + '.join(payment_equations)} = {total_revenue:.2f}"
    else:
        final_equation = f"Total Revenue ('{case_name}') = {total_revenue:.2f}"
    
    print(final_equation)
    print("-" * 30)
    
    return total_revenue

# Bids for the two scenarios
bids_x = [100, 20, 5]
bids_y = [100, 60, 10]

# Calculate x and y
x = calculate_auction_revenue(bids_x, 'x')
y = calculate_auction_revenue(bids_y, 'y')

# Print the final result tuple
print(f"\nThe final result for (x, y) is: ({int(x)}, {int(y)})")

<<< (10, 80) >>>