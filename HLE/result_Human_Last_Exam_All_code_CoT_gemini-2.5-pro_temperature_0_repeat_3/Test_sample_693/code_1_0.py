import collections

def solve_auction():
    """
    Solves the divisible item auction problem by calculating revenues x and y.
    """

    def get_allocation(bids_dict):
        """
        Calculates the allocation of the item based on the auction rule.
        The item is equally divided among bidders whose bid is at least 1/2 the highest bid.

        Args:
            bids_dict (dict): A dictionary of {bidder_id: bid_value}.

        Returns:
            dict: A dictionary of {bidder_id: allocation_share}.
        """
        if not bids_dict:
            return {}

        highest_bid = 0
        if bids_dict:
            highest_bid = max(bids_dict.values())

        threshold = highest_bid / 2.0
        winners = [bidder_id for bidder_id, bid in bids_dict.items() if bid >= threshold]

        allocation = collections.defaultdict(float)
        if winners:
            share = 1.0 / len(winners)
            for winner_id in winners:
                allocation[winner_id] = share
        
        return allocation

    def calculate_revenue(bids, case_name):
        """
        Calculates the total revenue for a given set of bids using the VCG payment rule.
        
        Args:
            bids (list): A list of bid values.
            case_name (str): The name of the case ('x' or 'y') for printing.

        Returns:
            float: The total revenue for the auction.
        """
        print(f"--- Calculating Revenue '{case_name}' for Bids: {bids} ---")
        
        # Use a dictionary to map bidder IDs to bids for clarity
        bids_dict = {i + 1: bid for i, bid in enumerate(bids)}
        
        # 1. Determine the outcome with all bidders present
        full_allocation = get_allocation(bids_dict)
        
        print(f"Initial allocation with all bidders:")
        highest_bid = max(bids_dict.values()) if bids_dict else 0
        winners = [bidder for bidder, share in full_allocation.items() if share > 0]
        print(f"  Highest Bid: {highest_bid}")
        print(f"  Winning Threshold ({highest_bid}/2): {highest_bid / 2.0}")
        if winners:
            print(f"  Winners: Bidder(s) {winners} (with bids {[bids_dict[w] for w in winners]})")
            print(f"  Allocation: Each winner gets a 1/{len(winners)} share.")
        else:
            print("  No winners.")
        
        total_revenue = 0
        payments = {}
        
        # 2. For each bidder, calculate their VCG payment
        for bidder_id in sorted(bids_dict.keys()):
            bidder_bid = bids_dict[bidder_id]
            
            # Create a dictionary of the other bidders' bids
            other_bidders_dict = bids_dict.copy()
            del other_bidders_dict[bidder_id]
            
            # a) Calculate the welfare of others WITH this bidder present
            welfare_others_with_bidder = sum(
                other_bidders_dict[other_id] * full_allocation[other_id]
                for other_id in other_bidders_dict
            )
            
            # b) Calculate the welfare of others WITHOUT this bidder
            alloc_without_bidder = get_allocation(other_bidders_dict)
            welfare_others_without_bidder = sum(
                other_bidders_dict[other_id] * alloc_without_bidder[other_id]
                for other_id in other_bidders_dict
            )
            
            # c) The payment is the "harm" caused to others
            payment = welfare_others_without_bidder - welfare_others_with_bidder
            payments[bidder_id] = payment
            total_revenue += payment
            
            print(f"\n  Calculating payment for Bidder {bidder_id} (bid={bidder_bid}):")
            print(f"    Welfare of others WITH Bidder {bidder_id}: {welfare_others_with_bidder:.2f}")
            print(f"    Welfare of others WITHOUT Bidder {bidder_id}: {welfare_others_without_bidder:.2f}")
            print(f"    Payment = {welfare_others_without_bidder:.2f} - {welfare_others_with_bidder:.2f} = {payment:.2f}")

        print("\nFinal Revenue Calculation:")
        payment_values = [payments[i] for i in sorted(payments.keys())]
        payment_strings = [f"{p:.2f}" for p in payment_values]
        print(f"  {case_name} = {' + '.join(payment_strings)} = {total_revenue:.2f}")
        return total_revenue

    # Case x: Bids are 100, 20, 5
    bids_x = [100, 20, 5]
    x = calculate_revenue(bids_x, 'x')
    
    print("\n" + "="*60 + "\n")
    
    # Case y: Bids are 100, 60, 10
    bids_y = [100, 60, 10]
    y = calculate_revenue(bids_y, 'y')
    
    print("\n" + "="*60 + "\n")
    print(f"The final result is (x, y) = ({int(x)}, {int(y)})")

solve_auction()