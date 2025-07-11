import math

def calculate_revenue(bids, case_label):
    """
    Calculates the total revenue for a given set of bids in a truthful, 
    divisible item auction, and prints the step-by-step calculation.

    The VCG payment mechanism is used to ensure truthfulness.
    """
    print(f"--- Calculating revenue '{case_label}' for bids: {bids} ---")
    
    n = len(bids)
    if n == 0:
        print(f"Revenue '{case_label}' = 0")
        return 0

    # Step 1: Determine winners and allocation with all bidders present
    b_max = max(bids)
    threshold = 0.5 * b_max
    
    winner_indices = [i for i, b in enumerate(bids) if b >= threshold]
    
    if not winner_indices:
        print("No bidders met the winning threshold.")
        print(f"Revenue '{case_label}' = 0")
        return 0
        
    k = len(winner_indices)
    allocation_per_winner = 1.0 / k

    total_revenue = 0.0
    payments = []
    payment_equations = []
    
    # Step 2: Calculate VCG payment for each winning bidder
    for i in winner_indices:
        # A) Calculate social welfare of other bidders *with* bidder i present
        welfare_others_with_i = 0.0
        for j in winner_indices:
            if i != j:
                welfare_others_with_i += bids[j] * allocation_per_winner

        # B) Calculate social welfare of other bidders if bidder i were absent
        bids_without_i = [bids[j] for j in range(n) if j != i]
        
        welfare_others_without_i = 0.0
        if bids_without_i:
            sub_b_max = max(bids_without_i)
            sub_threshold = 0.5 * sub_b_max
            
            sub_winners_bids = [b for b in bids_without_i if b >= sub_threshold]
            
            if sub_winners_bids:
                sub_k = len(sub_winners_bids)
                sub_alloc = 1.0 / sub_k
                welfare_others_without_i = sum(sub_winners_bids) * sub_alloc

        # C) The payment is the difference (the social cost imposed by bidder i)
        payment_i = welfare_others_without_i - welfare_others_with_i
        total_revenue += payment_i
        payments.append(payment_i)
        
        print(f"Bidder with bid {bids[i]} is a winner. Payment calculation:")
        print(f"  Payment = (Welfare of others without this bidder) - (Welfare of others with this bidder)")
        print(f"  Payment = {welfare_others_without_i:.2f} - {welfare_others_with_i:.2f} = {payment_i:.2f}")

    # Final revenue printout for the case
    payment_strings = [f"{p:.2f}" for p in payments]
    equation_str = " + ".join(payment_strings)
    
    print("\nTotal revenue equation:")
    if len(payments) > 1:
        print(f"Revenue '{case_label}' = {equation_str} = {total_revenue:.2f}")
    elif len(payments) == 1:
        print(f"Revenue '{case_label}' = {total_revenue:.2f}")
    else:
        print(f"Revenue '{case_label}' = 0.00")
        
    return total_revenue

if __name__ == "__main__":
    # Bids for case x
    bids_x = [100, 20, 5]
    x = calculate_revenue(bids_x, 'x')

    print("\n" + "="*50 + "\n")

    # Bids for case y
    bids_y = [100, 60, 10]
    y = calculate_revenue(bids_y, 'y')

    # Final combined answer
    print("\n" + "="*50 + "\n")
    print("### Final Answer ###")
    print(f"The revenue x is {x:.2f}")
    print(f"The revenue y is {y:.2f}")
    print(f"The result (x, y) is: ({x:.2f}, {y:.2f})")
    # The final tuple format for the platform.
    # <<<(20.0, 80.0)>>>