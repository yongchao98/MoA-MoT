import math

def calculate_welfare_and_winners(bids):
    """
    Helper function to calculate winners and social welfare for a given set of bids.
    Social Welfare = sum of (1/k) * bid for each of k winners.
    """
    if not bids:
        return 0, [], 0

    highest_bid = max(bids)
    threshold = 0.5 * highest_bid
    
    winners = [b for b in bids if b >= threshold]
    k = len(winners)

    if k == 0:
        return 0, [], 0

    # For a divisible good, social welfare is the sum of valuations of the winners' shares
    welfare = (1/k) * sum(winners)
    return welfare, winners, k

def calculate_vcg_revenue(bids, case_name):
    """
    Calculates the total revenue of the auction using VCG payments.
    """
    print(f"--- Calculating revenue {case_name} for bids {bids} ---")
    
    # Step 1: Analyze the outcome with all bidders
    original_welfare, original_winners, k = calculate_welfare_and_winners(bids)
    
    if k == 0:
        print("There are no winners, so the revenue is 0.")
        return 0

    highest_bid = max(bids)
    threshold = 0.5 * highest_bid
    print(f"1. With bids {bids}:")
    print(f"   - The highest bid is {highest_bid:.0f}.")
    print(f"   - The winning threshold is 0.5 * {highest_bid:.0f} = {threshold:.0f}.")
    print(f"   - Winners (bids >= {threshold:.0f}) are: {original_winners}")
    print(f"   - The item is split among {k} winner(s).\n")

    total_revenue = 0
    payment_components = []
    
    # Step 2: Calculate payment for each winner
    print("2. Calculating VCG payments for each winner:")
    for i, winner_bid in enumerate(original_winners):
        # Create the list of bids without the current winner
        bids_without_winner = list(bids)
        bids_without_winner.remove(winner_bid)
        
        # Calculate what the welfare would have been without this winner
        welfare_without, _, _ = calculate_welfare_and_winners(bids_without_winner)
        
        # Calculate the welfare of OTHERS in the original auction
        welfare_of_others_in_original_auction = original_welfare - (1/k * winner_bid)
        
        # VCG Payment = (Welfare without me) - (Welfare of others with me)
        payment = welfare_without - welfare_of_others_in_original_auction
        payment_components.append(payment)
        total_revenue += payment
        
        print(f"   - For winner with bid {winner_bid:.0f}:")
        print(f"     a) If this bidder did not participate, the bids would be {bids_without_winner}.")
        highest_bid_without = max(bids_without_winner) if bids_without_winner else 0
        threshold_without = 0.5 * highest_bid_without
        winners_without = [b for b in bids_without_winner if b >= threshold_without]
        print(f"        The winners would be {winners_without}, and welfare would be {welfare_without:.0f}.")
        print(f"     b) In the original auction, the welfare of OTHERS was {welfare_of_others_in_original_auction:.0f}.")
        print(f"     c) Payment = (welfare without them) - (welfare of others with them)")
        print(f"        Payment = {welfare_without:.0f} - {welfare_of_others_in_original_auction:.0f} = {payment:.0f}.\n")

    # Step 3: Print total revenue
    print(f"3. Total revenue {case_name} is the sum of payments.")
    payment_str = " + ".join([f"{p:.0f}" for p in payment_components])
    print(f"   Total revenue {case_name} = {payment_str} = {total_revenue:.0f}.\n")
    return total_revenue

def solve():
    """
    Main function to solve the problem for x and y.
    """
    bids_x = [100, 20, 5]
    bids_y = [100, 60, 10]

    x = calculate_vcg_revenue(bids_x, 'x')
    y = calculate_vcg_revenue(bids_y, 'y')

    print("-----------------------------------------")
    print("Final Result:")
    print(f"The revenue x is {x:.0f}.")
    print(f"The revenue y is {y:.0f}.")
    print(f"The final answer for (x, y) is ({x:.0f}, {y:.0f}).")


solve()
<<< (20, 80) >>>