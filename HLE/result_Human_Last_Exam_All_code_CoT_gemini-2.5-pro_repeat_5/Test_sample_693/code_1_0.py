import collections

def calculate_vcg_revenue(bids, case_name):
    """
    Calculates the total revenue of the specified divisible good auction
    using the VCG payment rule.

    The function also prints the step-by-step calculation.
    """
    print(f"--- Calculating revenue '{case_name}' for bids: {bids} ---")

    if not bids:
        print("No bids, revenue is 0.")
        return 0

    # Step 1: Determine the outcome with all bidders present.
    # We use a dictionary to handle potential non-unique bids and track bidders.
    bidders = collections.OrderedDict(enumerate(bids))
    
    highest_bid = float(max(bidders.values()))
    threshold = highest_bid / 2.0
    
    # Use original indices to identify winners
    winner_indices = {i for i, bid in bidders.items() if bid >= threshold}
    k = len(winner_indices)

    if k == 0:
        print(f"Highest bid is {highest_bid}, threshold is {threshold}.")
        print("There are no winners.")
        print(f"Total revenue {case_name} = 0.0")
        return 0

    winner_bids = [bidders[i] for i in winner_indices]
    print(f"Highest bid is {highest_bid}, threshold is {threshold}.")
    print(f"Winners are bidders with bids {set(winner_bids)} ({k} winners). Each gets a 1/{k} share.")
    print("-" * 20)

    total_revenue = 0
    winner_payments = []

    # Step 2: For each winner, calculate their VCG payment.
    for i in winner_indices:
        winner_bid = bidders[i]
        
        # Calculate the social welfare of OTHERS, with winner `i` present.
        welfare_of_others_with_i = 0
        for j in winner_indices:
            if i == j:
                continue
            welfare_of_others_with_i += bidders[j] * (1.0 / k)

        # Calculate the social welfare in a world WITHOUT winner `i`.
        bidders_without_i = {idx: b for idx, b in bidders.items() if idx != i}
        
        welfare_without_i = 0
        if bidders_without_i:
            highest_bid_without_i = float(max(bidders_without_i.values()))
            threshold_without_i = highest_bid_without_i / 2.0
            
            winners_without_i_indices = {idx for idx, b in bidders_without_i.items() if b >= threshold_without_i}
            k_without_i = len(winners_without_i_indices)

            if k_without_i > 0:
                for j in winners_without_i_indices:
                    welfare_without_i += bidders[j] * (1.0 / k_without_i)
        
        # VCG Payment = (Welfare without i) - (Welfare of others with i)
        payment = welfare_without_i - welfare_of_others_with_i
        winner_payments.append(payment)
        total_revenue += payment
        
        # Print the detailed calculation for this winner's payment
        print(f"Payment for bidder with bid {winner_bid}:")
        print(f"  Welfare of others (with this bidder) = {welfare_of_others_with_i:.1f}")
        print(f"  Welfare of world (without this bidder) = {welfare_without_i:.1f}")
        print(f"  Payment = {welfare_without_i:.1f} - {welfare_of_others_with_i:.1f} = {payment:.1f}")
        print("-" * 20)

    payment_summary = " + ".join([f"{p:.1f}" for p in winner_payments])
    print(f"Total revenue {case_name} = {payment_summary} = {total_revenue:.1f}\n")
    return total_revenue

if __name__ == "__main__":
    bids_x = [100, 20, 5]
    x = calculate_vcg_revenue(bids_x, 'x')

    bids_y = [100, 60, 10]
    y = calculate_vcg_revenue(bids_y, 'y')

    print("--- Final Answer ---")
    # The final result is a tuple (x, y)
    final_answer = (x, y)
    print(f"The result (x, y) is: {final_answer}")
    print(f"<<<{final_answer}>>>")
