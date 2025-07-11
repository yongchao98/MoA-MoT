import math

def calculate_vcg_revenue(bids, label):
    """
    Calculates the total revenue of a divisible item auction using VCG payments.

    Args:
        bids (list): A list of the bids from all bidders.
        label (str): The label for the revenue ('x' or 'y') for printing.
    """
    print(f"--- Calculating revenue {label} for bids {tuple(bids)} ---")

    if not bids:
        print("No bids provided. Revenue is 0.")
        return 0

    # Find the set of winners with all bidders present
    b_max = max(bids) if bids else 0
    if b_max <= 0:
        print("Highest bid is not positive. No winners.")
        return 0
    
    threshold = 0.5 * b_max
    # Keep track of original indices
    indexed_bids = list(enumerate(bids))
    
    winners = [b for b in indexed_bids if b[1] >= threshold]
    k = len(winners)

    if k == 0:
        print(f"Highest bid is {b_max}. Threshold is {threshold}. There are no winners.")
        print(f"Revenue {label} = 0")
        return 0
    
    winning_bids = [w[1] for w in winners]
    print(f"Highest bid is {b_max}. Winning threshold is {b_max} * 0.5 = {threshold}.")
    print(f"The winners are bidders with bids {winning_bids}, who get 1/{k} of the item each.")
    print("-" * 20)

    total_revenue = 0
    winner_payments = []

    # Calculate payment for each winner
    for winner_idx, winner_bid in winners:
        
        # Calculate social value of OTHERS when this winner is PRESENT
        value_of_others_present = 0
        for other_winner_idx, other_winner_bid in winners:
            if other_winner_idx != winner_idx:
                value_of_others_present += other_winner_bid / k

        # Calculate social value if this winner were ABSENT
        bids_without_winner = [b[1] for b in indexed_bids if b[0] != winner_idx]
        value_without_winner = 0
        
        # Determine winners and social value in the world without this bidder
        sw_details = ""
        if bids_without_winner:
            b_max_prime = max(bids_without_winner)
            threshold_prime = 0.5 * b_max_prime
            winners_prime = [b for b in bids_without_winner if b >= threshold_prime]
            k_prime = len(winners_prime)
            if k_prime > 0:
                value_without_winner = sum(winners_prime) / k_prime
            sw_details = (f"If this bidder were absent, bids would be {tuple(bids_without_winner)}. "
                          f"Highest bid would be {b_max_prime}, threshold {threshold_prime}. "
                          f"The winners would be {winners_prime}, so social value = {value_without_winner:.2f}.")
        else:
            sw_details = "If this bidder were absent, there would be no other bidders, so social value = 0."
            
        
        payment = value_without_winner - value_of_others_present
        total_revenue += payment
        winner_payments.append(str(int(round(payment))))
        
        print(f"Payment for bidder with bid {winner_bid}:")
        print(f"- {sw_details}")
        print(f"- When this bidder is present, the value allocated to others is {value_of_others_present:.2f}.")
        print(f"- Payment = {value_without_winner:.2f} - {value_of_others_present:.2f} = {payment:.2f}")
        print()

    print(f"Total revenue {label} = {' + '.join(winner_payments)} = {int(round(total_revenue))}")
    return total_revenue

if __name__ == '__main__':
    bids_x = [100, 20, 5]
    x = calculate_vcg_revenue(bids_x, 'x')

    print("\n============================================\n")

    bids_y = [100, 60, 10]
    y = calculate_vcg_revenue(bids_y, 'y')

    print("\n============================================\n")
    print(f"The final answer is (x, y) = ({int(round(x))}, {int(round(y))})")
    print("<<<({}, {})>>>".format(int(round(x)), int(round(y))))
