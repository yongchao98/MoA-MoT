def calculate_auction_revenue(bids, case_name):
    """
    Calculates and prints the details of an auction's revenue.

    Args:
        bids (list): A list of numbers representing the bids.
        case_name (str): The name of the variable for the revenue (e.g., 'x' or 'y').
    """
    print(f"--- Calculating revenue {case_name} for bids: {bids} ---")

    if not bids:
        print("No bids provided. Revenue is 0.")
        return 0

    h = max(bids)
    threshold = h / 2.0
    print(f"The highest bid is {h}, so the winning threshold is {h}/2 = {threshold}.")

    winners = []
    # Use enumerate to keep track of the original index of each bid
    for i, bid in enumerate(bids):
        if bid >= threshold:
            winners.append({'id': i, 'bid': bid})

    if not winners:
        print("There are no winners in this auction.")
        print(f"Total revenue {case_name} = 0.0")
        return 0

    num_winners = len(winners)
    winner_bids_str = ", ".join([str(w['bid']) for w in winners])
    print(f"The winners are the bidders with bids: {winner_bids_str}. The number of winners is {num_winners}.")

    total_revenue = 0
    payment_expressions = []

    for winner in winners:
        # Get all bids except the current winner's bid
        other_bids = [b for i, b in enumerate(bids) if i != winner['id']]
        
        # h_minus_i is the highest bid among the others
        h_minus_i = max(other_bids) if other_bids else 0
        
        # The price per unit for this winner
        price_per_unit = h_minus_i / 2.0
        
        # The share of the item this winner receives
        share = 1.0 / num_winners
        
        # The payment for this winner
        payment = share * price_per_unit
        
        print(f"Payment for the winner with bid {winner['bid']}:")
        # Print the equation with all the numbers
        print(f"  = (1/{num_winners}) * (max({other_bids}) / 2)")
        print(f"  = (1/{num_winners}) * ({h_minus_i} / 2)")
        print(f"  = {share:.2f} * {price_per_unit:.2f} = {payment:.2f}")

        total_revenue += payment
        payment_expressions.append(f"{payment:.2f}")

    revenue_equation = " + ".join(payment_expressions)
    print(f"\nThe total revenue {case_name} = {revenue_equation} = {total_revenue:.2f}")
    return total_revenue

if __name__ == "__main__":
    # Case for x: bids are 100, 20, and 5
    bids_x = [100, 20, 5]
    x = calculate_auction_revenue(bids_x, 'x')
    
    print("\n" + "="*50 + "\n")
    
    # Case for y: bids are 100, 60, and 10
    bids_y = [100, 60, 10]
    y = calculate_auction_revenue(bids_y, 'y')
    
    print("\n" + "="*50 + "\n")
    
    # Final Answer
    print(f"The final result is (x, y) = ({x:.2f}, {y:.2f})")