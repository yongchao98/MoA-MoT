def calculate_auction_outcome(bids):
    """
    Determines winners and their allocations based on the auction rules.
    Returns a list of winner indices and a corresponding list of allocations.
    """
    if not bids:
        return [], []
    
    max_bid = 0
    # Handle the case where all bids are 0 or negative
    if any(b > 0 for b in bids):
        max_bid = max(bids)

    threshold = max_bid / 2.0
    
    winner_indices = [i for i, b in enumerate(bids) if b >= threshold]
    num_winners = len(winner_indices)
    
    allocations = [0.0] * len(bids)
    if num_winners > 0:
        share = 1.0 / num_winners
        for i in winner_indices:
            allocations[i] = share
            
    return winner_indices, allocations

def calculate_revenue(bids):
    """
    Calculates the total revenue from an auction using VCG payments.
    """
    # 1. Calculate the outcome with all bidders present
    original_winner_indices, original_allocations = calculate_auction_outcome(bids)
    
    if not original_winner_indices:
        return 0, []
        
    # Calculate the utility (value) each bidder gets in the original auction
    original_utilities = [b * q for b, q in zip(bids, original_allocations)]
    original_total_utility = sum(original_utilities)
    
    total_revenue = 0
    payments = []
    
    # 2. For each winner, calculate their VCG payment
    for i in original_winner_indices:
        # This is the value all other bidders get when bidder `i` participates
        utility_of_others_with_i = original_total_utility - original_utilities[i]
        
        # 3. Create a hypothetical auction without bidder `i`
        bids_without_i = [b for j, b in enumerate(bids) if j != i]
        
        # Calculate the outcome and utility for this hypothetical auction
        _, hypothetical_allocations = calculate_auction_outcome(bids_without_i)
        
        # This is the value the other bidders would have gotten without bidder `i`
        utility_of_others_without_i = sum(b * q for b, q in zip(bids_without_i, hypothetical_allocations))
        
        # 4. The payment is the difference ("harm" caused by i's participation)
        payment_i = utility_of_others_without_i - utility_of_others_with_i
        payments.append(payment_i)
        total_revenue += payment_i
        
    return total_revenue, payments

def main():
    """
    Solves the user's request for auction revenues (x, y).
    """
    bids_x = [100, 20, 5]
    bids_y = [100, 60, 10]

    # Calculate x
    x, payments_x = calculate_revenue(bids_x)
    print(f"For bids {bids_x}:")
    payment_str_x = " + ".join([f"{p:.0f}" for p in payments_x])
    print(f"Revenue x from payments: {payment_str_x} = {x:.0f}")

    print("-" * 20)

    # Calculate y
    y, payments_y = calculate_revenue(bids_y)
    print(f"For bids {bids_y}:")
    payment_str_y = " + ".join([f"{p:.0f}" for p in payments_y])
    print(f"Revenue y from payments: {payment_str_y} = {y:.0f}")
    
    print("-" * 20)
    
    # Final answer
    print(f"The final answer for (x, y) is ({x:.0f}, {y:.0f})")

if __name__ == "__main__":
    main()