def calculate_revenue(bids, case_name):
    """
    Calculates the revenue of a divisible item auction using VCG payments.
    It also prints the step-by-step calculation.
    """
    print(f"--- Calculating revenue '{case_name}' for bids: {bids} ---")
    
    # Use a dictionary to keep track of original bidder indices
    bidders = {i: bid for i, bid in enumerate(bids)}

    # Step 1: Identify winners with all bidders present
    if not bidders:
        print("No bidders, revenue is 0.")
        return 0
        
    max_bid = max(bidders.values())
    threshold = 0.5 * max_bid
    
    winners = {i: bid for i, bid in bidders.items() if bid >= threshold}
    num_winners = len(winners)

    print(f"The highest bid is {max_bid}, so the winning threshold is 0.5 * {max_bid} = {threshold}.")
    
    if not winners:
        print("There are no winners. Total revenue is 0.\n")
        return 0
    
    winner_bids = list(winners.values())
    print(f"The winners are the bidders with bids {winner_bids}.")
    print(f"The item is split {num_winners} way(s).\n")

    total_revenue = 0
    
    # Step 2 & 3: Calculate VCG payment for each winner and sum them up
    for winner_idx, winner_bid in winners.items():
        print(f"Calculating payment for the winner with bid {winner_bid}:")
        
        # Calculate social welfare of OTHERS with the winner present
        welfare_of_others_with_winner = 0
        for other_idx, other_bid in winners.items():
            if other_idx != winner_idx:
                welfare_of_others_with_winner += other_bid * (1 / num_winners)
        
        print(f"  - Social welfare of other winners (if any) with this winner present:")
        if num_winners > 1:
            others_explanation = " + ".join([f"{bid} * 1/{num_winners}" for bid in winners.values() if bid != winner_bid])
            print(f"    = {others_explanation} = {welfare_of_others_with_winner:.2f}")
        else:
            print("    = 0 (since there are no other winners)")

        # Calculate social welfare WITHOUT the winner
        bids_without_winner = bidders.copy()
        del bids_without_winner[winner_idx]
        
        welfare_without_winner = 0
        if bids_without_winner:
            new_max_bid = max(bids_without_winner.values())
            new_threshold = 0.5 * new_max_bid
            new_winners = {i: bid for i, bid in bids_without_winner.items() if bid >= new_threshold}
            num_new_winners = len(new_winners)
            
            if new_winners:
                for bid in new_winners.values():
                    welfare_without_winner += bid * (1 / num_new_winners)
                
                print(f"  - Social welfare if this winner had not participated:")
                print(f"    The remaining bids would be {list(bids_without_winner.values())}.")
                print(f"    The new highest bid would be {new_max_bid}, and the threshold {new_threshold}.")
                welfare_explanation = " + ".join([f"{bid} * 1/{num_new_winners}" for bid in new_winners.values()])
                print(f"    The new welfare would be: {welfare_explanation} = {welfare_without_winner:.2f}")
            else:
                print("  - If this winner had not participated, there would be no winners, so social welfare is 0.")
        else:
            print("  - If this winner had not participated, there would be no bidders, so social welfare is 0.")
            
        # Calculate the payment
        payment = welfare_without_winner - welfare_of_others_with_winner
        total_revenue += payment
        
        print(f"  - The payment is calculated as: (welfare without winner) - (welfare of others with winner)")
        print(f"    Payment = {welfare_without_winner:.2f} - {welfare_of_others_with_winner:.2f} = {payment:.2f}\n")

    print(f"Total revenue for case '{case_name}' is the sum of all payments: {total_revenue:.2f}\n")
    return total_revenue

if __name__ == "__main__":
    # Case x
    bids_x = [100, 20, 5]
    x = calculate_revenue(bids_x, 'x')
    
    # Case y
    bids_y = [100, 60, 10]
    y = calculate_revenue(bids_y, 'y')
    
    # Final Answer
    print("--- Final Answer ---")
    print(f"The revenue for x is {int(x)}.")
    print(f"The revenue for y is {int(y)}.")
    print(f"The result (x, y) is: ({int(x)}, {int(y)})")
    
    final_answer = (int(x), int(y))
    # The final answer is wrapped in <<<>>> as requested.
    # print(f"<<<{final_answer}>>>")