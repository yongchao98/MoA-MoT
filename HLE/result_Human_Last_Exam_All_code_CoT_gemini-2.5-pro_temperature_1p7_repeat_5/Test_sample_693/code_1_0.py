def calculate_revenue(bids: list) -> float:
    """
    Calculates the total revenue of the auction based on the VCG payment rule.
    """
    if not bids:
        return 0.0

    # Determine original winners
    b_max = max(bids)
    threshold = 0.5 * b_max
    # Assuming unique bids as in the problem, this identifies each winner
    winners = sorted([b for b in bids if b >= threshold], reverse=True)
    
    k = len(winners)
    if k == 0:
        return 0.0

    total_revenue = 0.0
    for winner_bid in winners:
        # Calculate the welfare of other participants WITH this winner
        # This is the sum of valuations of other winners for their shares
        welfare_of_others_with = (sum(winners) - winner_bid) / k
        
        # Create a list of bids for the hypothetical auction without this winner
        bids_without = [b for b in bids if b != winner_bid]
        
        # Calculate the total welfare of participants in this hypothetical auction
        hypothetical_welfare = 0.0
        if bids_without:
            b_max_prime = max(bids_without)
            threshold_prime = 0.5 * b_max_prime
            winners_prime = [b for b in bids_without if b >= threshold_prime]
            k_prime = len(winners_prime)
            if k_prime > 0:
                hypothetical_welfare = sum(winners_prime) / k_prime

        # The VCG payment is the harm caused to others
        payment = hypothetical_welfare - welfare_of_others_with
        total_revenue += payment
        
    return total_revenue

def solve():
    """
    Solves the problem by calculating revenues x and y for the two scenarios.
    """
    # Case x: Bids are 100, 20, and 5
    bids_x = [100, 20, 5]
    x = calculate_revenue(bids_x)

    # Case y: Bids are 100, 60, and 10
    bids_y = [100, 60, 10]
    y = calculate_revenue(bids_y)

    # The problem asks to output each number in the final equation.
    # We will print the value of x and y separately and then the final tuple (x, y).
    print(f"The revenue x for bids {bids_x} is {int(x)}")
    print(f"The revenue y for bids {bids_y} is {int(y)}")
    print(f"The final answer for (x, y) is ({int(x)}, {int(y)})")

solve()