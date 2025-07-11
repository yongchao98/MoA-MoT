import math

def calculate_outcome(bids_dict):
    """
    Calculates the winners, their shares, and the total welfare for a given set of bids.
    bids_dict is a dictionary of {bidder_id: bid_value}.
    """
    if not bids_dict:
        return {}, {}, 0.0

    highest_bid = 0
    for bid in bids_dict.values():
        if bid > highest_bid:
            highest_bid = bid

    threshold = highest_bid / 2.0
    
    winners = {bidder_id: bid for bidder_id, bid in bids_dict.items() if bid >= threshold}
    
    num_winners = len(winners)
    if num_winners == 0:
        return {}, {}, 0.0

    share_per_winner = 1.0 / num_winners
    
    shares = {bidder_id: share_per_winner for bidder_id in winners}
    
    total_welfare = sum(bid * share_per_winner for bid in winners.values())
    
    return winners, shares, total_welfare

def calculate_revenue_and_print(bids_list, case_name):
    """
    Calculates the total revenue for a list of bids based on the VCG mechanism
    and prints the detailed steps.
    """
    print(f"--- Calculating revenue '{case_name}' for bids: {bids_list} ---")
    
    # Create a dictionary with unique IDs for each bidder
    bids_dict = {i + 1: bid for i, bid in enumerate(bids_list)}

    main_winners, main_shares, main_welfare = calculate_outcome(bids_dict)
    
    winning_bids = sorted(list(main_winners.values()), reverse=True)
    print(f"Highest Bid: {max(bids_list) if bids_list else 0}")
    print(f"Winning Threshold (1/2 * Highest Bid): {(max(bids_list) if bids_list else 0) / 2.0}")
    print(f"Winning bids: {winning_bids}")
    if not winning_bids:
        print(f"Total revenue {case_name} = 0\n")
        return 0

    total_revenue = 0
    payments = []
    
    # Sort by bidder ID to ensure consistent output
    sorted_winner_ids = sorted(main_winners.keys())

    for bidder_id in sorted_winner_ids:
        bid_value = main_winners[bidder_id]
        share = main_shares[bidder_id]
        
        # Calculate welfare of others WITH this bidder
        welfare_of_others_with_bidder = main_welfare - (bid_value * share)

        # Calculate outcome WITHOUT this bidder
        bids_without_bidder = bids_dict.copy()
        del bids_without_bidder[bidder_id]
        _, _, welfare_without_bidder = calculate_outcome(bids_without_bidder)
        
        # VCG Payment
        payment = welfare_without_bidder - welfare_of_others_with_bidder
        total_revenue += payment
        payments.append(payment)
        
        print(f"Payment for bidder with bid {bid_value}: {welfare_without_bidder:.2f} - {welfare_of_others_with_bidder:.2f} = {payment:.2f}")

    payment_equation = " + ".join([f"{p:.2f}" for p in payments])
    print(f"Total revenue {case_name} = {payment_equation} = {total_revenue:.2f}\n")
    return total_revenue

if __name__ == "__main__":
    bids_x = [100, 20, 5]
    x = calculate_revenue_and_print(bids_x, 'x')
    
    bids_y = [100, 60, 10]
    y = calculate_revenue_and_print(bids_y, 'y')
    
    print("--- Final Answer ---")
    # Use math.floor to ensure integer output for clean formatting if they are whole numbers
    final_tuple = (int(math.floor(x)), int(math.floor(y)))
    print(f"(x, y) = {final_tuple}")
    print(f"<<<{final_tuple}>>>")
