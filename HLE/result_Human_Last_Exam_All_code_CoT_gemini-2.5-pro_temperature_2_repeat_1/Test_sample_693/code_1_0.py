import math

def get_auction_outcome(bids):
    """
    Calculates winners, allocation, and social welfare for a set of bids.
    """
    if not bids:
        return {}, 0, 0
    
    b_max = max(bids)
    threshold = 0.5 * b_max
    
    # Using a dictionary to handle potential duplicate bids while keeping them distinct
    winners = {i: bid for i, bid in enumerate(bids) if bid >= threshold}
    
    k = len(winners)
    if k == 0:
        return {}, 0, 0
    
    sw = sum(bid * (1/k) for bid in winners.values())
    return winners, k, sw

def calculate_revenue_and_print(bids, name):
    """
    Calculates the revenue and prints the breakdown.
    """
    print(f"--- Calculating revenue {name} for bids {bids} ---")
    
    # Original auction
    winners, k, sw_total = get_auction_outcome(bids)
    
    if k == 0:
        print(f"No winners. Revenue {name} = 0")
        return 0

    payments = []
    winner_bids = list(winners.values())
    print(f"Winning bids: {winner_bids}, each gets 1/{k} of the item.")
    
    # Calculate payment for each winner
    for winner_id, winner_bid in winners.items():
        # SW of others with this winner participating
        sw_others = sw_total - winner_bid * (1/k)
        
        # SW if this winner did not participate
        bids_without_winner = [b for i, b in enumerate(bids) if i != winner_id]
        _, _, sw_without_winner = get_auction_outcome(bids_without_winner)
        
        payment = sw_without_winner - sw_others
        payments.append(payment)
    
    revenue = sum(payments)
    
    # Format payment equation string
    payment_str = " + ".join([f"{p:.0f}" for p in payments])
    print(f"Revenue equation: {name} = {payment_str}")
    print(f"Total revenue {name} = {revenue:.0f}\n")
    
    return revenue

# Calculate x
bids_x = [100, 20, 5]
x = calculate_revenue_and_print(bids_x, 'x')

# Calculate y
bids_y = [100, 60, 10]
y = calculate_revenue_and_print(bids_y, 'y')

# Final answer
print("--------------------")
print(f"The final result (x, y) is: ({x:.0f}, {y:.0f})")
