import math

def calculate_revenue(bids, case_name):
    """
    Calculates the total auction revenue based on the specified allocation rule
    and a truthful (VCG) payment mechanism.

    Args:
        bids (list): A list of the bids from all bidders.
        case_name (str): The name of the case ('x' or 'y') for printing.
    """
    print(f"--- Calculating revenue for case '{case_name}' with bids {bids} ---")

    if not bids:
        print(f"Total revenue for '{case_name}' is: 0\n")
        return 0

    # 1. Determine the outcome with all bidders present (the "main auction")
    highest_bid = max(bids)
    threshold = highest_bid / 2.0
    
    main_winners = []
    for bid in bids:
        if bid >= threshold:
            main_winners.append(bid)
    num_main_winners = len(main_winners)

    if num_main_winners == 0:
        print("There are no winners in the auction.")
        print(f"Total revenue for '{case_name}' is: 0\n")
        return 0

    print(f"In the main auction, the highest bid is {highest_bid}, threshold is {threshold}.")
    print(f"The winners are bidders with bids {main_winners}, and the item is split {num_main_winners} ways.")

    total_revenue = 0

    # 2. Iterate through each bidder to calculate their VCG payment
    for i in range(len(bids)):
        bidder_bid = bids[i]

        # Losers pay nothing. We only need to calculate payments for winners.
        if bidder_bid not in main_winners:
            continue

        # 3. Calculate the total welfare of OTHERS with bidder i present
        welfare_of_others_with_i = 0
        for j in range(len(bids)):
            if i == j:
                continue
            # Check if bidder j is a winner in the main auction
            if bids[j] in main_winners:
                welfare_of_others_with_i += bids[j] / num_main_winners

        # 4. Calculate the total welfare of OTHERS with bidder i absent
        bids_without_i = bids[:i] + bids[i+1:]
        
        welfare_of_others_without_i = 0
        if bids_without_i:
            highest_bid_without_i = max(bids_without_i)
            threshold_without_i = highest_bid_without_i / 2.0
            
            winners_without_i = []
            for bid in bids_without_i:
                if bid >= threshold_without_i:
                    winners_without_i.append(bid)
            
            num_winners_without_i = len(winners_without_i)

            if num_winners_without_i > 0:
                # The total value for the winners is the sum of their bids (valuations)
                # divided by the number of winners (since the item is split equally).
                welfare_of_others_without_i = sum(winners_without_i) / num_winners_without_i

        # 5. The payment is the "harm" caused to others
        payment_i = welfare_of_others_without_i - welfare_of_others_with_i
        
        # The VCG payment should not be negative.
        payment_i = max(0, payment_i)

        print(f"Payment for winner with bid {bidder_bid}: (Welfare of others without them) - (Welfare of others with them)")
        print(f"  = {welfare_of_others_without_i:.2f} - {welfare_of_others_with_i:.2f} = {payment_i:.2f}")
        
        total_revenue += payment_i

    print(f"Total revenue for '{case_name}' is: {total_revenue:.2f}\n")
    return total_revenue

if __name__ == "__main__":
    # Bids for the two scenarios
    bids_x = [100, 20, 5]
    bids_y = [100, 60, 10]

    # Calculate revenues
    x = calculate_revenue(bids_x, 'x')
    y = calculate_revenue(bids_y, 'y')

    # Print the final result
    # The problem asks for integers, so we format them as such.
    final_x = int(round(x))
    final_y = int(round(y))
    
    print("--- Final Result ---")
    print(f"The revenue x is {final_x}.")
    print(f"The revenue y is {final_y}.")
    print(f"The final answer for (x, y) is ({final_x}, {final_y}).")
    print(f"<<<({final_x}, {final_y})>>>")