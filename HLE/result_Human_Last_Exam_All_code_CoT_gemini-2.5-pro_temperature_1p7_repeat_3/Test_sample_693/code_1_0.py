def calculate_auction_revenue(bids, scenario_name):
    """
    Calculates the revenue of a single-item divisible auction based on the given rules.

    Args:
        bids (list): A list of numerical bids.
        scenario_name (str): The name of the scenario ('x' or 'y') for printing.

    Returns:
        float: The total revenue from the auction.
    """
    print(f"--- Calculating revenue {scenario_name} for bids: {bids} ---")
    
    if not bids:
        print("No bids provided. Revenue is 0.")
        return 0

    b_max = max(bids)
    winning_threshold = 0.5 * b_max
    
    print(f"Highest bid (b_max): {b_max}")
    print(f"Winning threshold (0.5 * b_max): {winning_threshold}")
    print(f"Winners are bidders with bids >= {winning_threshold}")

    winners = []
    for i, bid in enumerate(bids):
        if bid >= winning_threshold:
            winners.append((i, bid))

    total_revenue = 0
    payment_components = []

    if not winners:
        print("There are no winners.")
    else:
        for i, bid in winners:
            # Create a list of other bids
            other_bids = bids[:i] + bids[i+1:]
            
            # The payment for a winner is based on the highest bid among others
            b_max_other = max(other_bids) if other_bids else 0
            payment = 0.5 * b_max_other
            
            print(f"Winner with bid {bid} pays: 0.5 * max(other bids {other_bids}) = {payment}")
            
            total_revenue += payment
            payment_components.append(str(payment))
    
    # Show the final equation for the revenue
    if len(payment_components) > 1:
        equation = " + ".join(payment_components)
        print(f"Total revenue {scenario_name} = {equation} = {total_revenue}")
    else:
        print(f"Total revenue {scenario_name} = {total_revenue}")
        
    print("-" * (len(f"--- Calculating revenue {scenario_name} for bids: {bids} ---")))
    return total_revenue

if __name__ == "__main__":
    bids_x = [100, 20, 5]
    bids_y = [100, 60, 10]
    
    x = calculate_auction_revenue(bids_x, 'x')
    print()
    y = calculate_auction_revenue(bids_y, 'y')
    
    print(f"\nFinal result: (x, y) = ({x}, {y})")
    
    # Final answer in the specified format
    # The output format is <<<answer content>>>. The content is the tuple (x, y).
    # Since the values are integers, we use (10, 80).
    final_answer = (int(x), int(y))
    # print(f"<<<{final_answer}>>>") # This would be the submission format line