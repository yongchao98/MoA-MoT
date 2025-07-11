import collections

def get_auction_outcome(bids):
    """
    Determines winners and welfare for a given set of bids.
    """
    if not bids:
        return [], 0

    max_bid = max(bids)
    threshold = 0.5 * max_bid
    
    winners = [b for b in bids if b >= threshold]
    
    if not winners:
        return [], 0
        
    num_winners = len(winners)
    social_welfare = sum(winners) / num_winners
    
    return winners, social_welfare

def calculate_revenue(bids, case_name):
    """
    Calculates the total revenue of the auction using VCG payments.
    """
    print(f"--- Calculating Revenue '{case_name}' for bids {bids} ---")
    
    # Original auction outcome
    original_winners, _ = get_auction_outcome(bids)
    
    if not original_winners:
        print("No winners in the auction.")
        print(f"Total Revenue '{case_name}': 0")
        return 0
        
    print(f"Original winners are: {original_winners}")
    num_original_winners = len(original_winners)
    
    total_revenue = 0
    
    # Use Counter to handle duplicate bids correctly
    original_bids_counts = collections.Counter(bids)

    # Calculate payment for each winner
    for winner_bid in original_winners:
        # Create the list of bids without one instance of the current winner's bid
        bids_without_winner = list(bids)
        bids_without_winner.remove(winner_bid)
        
        # Calculate welfare if this winner had not participated
        _, welfare_without_winner = get_auction_outcome(bids_without_winner)
        
        # Calculate welfare of OTHERS in the original auction
        other_winners = list(original_winners)
        other_winners.remove(winner_bid)
        
        welfare_of_others = 0
        if other_winners:
            welfare_of_others = sum(other_winners) / num_original_winners
            
        payment = welfare_without_winner - welfare_of_others
        
        print(f"Payment for bidder with bid {winner_bid}:")
        print(f"  Welfare without this bidder = {welfare_without_winner}")
        print(f"  Welfare of other winners in original auction = {welfare_of_others}")
        print(f"  Payment = {welfare_without_winner} - {welfare_of_others} = {payment}")
        
        total_revenue += payment
        
    print(f"Total Revenue '{case_name}': {total_revenue}\n")
    return total_revenue

def main():
    bids1 = [100, 20, 5]
    bids2 = [100, 60, 10]

    x = calculate_revenue(bids1, 'x')
    y = calculate_revenue(bids2, 'y')

    print("--- Final Result ---")
    print(f"The value of (x, y) is: ({x}, {y})")

if __name__ == "__main__":
    main()