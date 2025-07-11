def solve_divisible_item_auction():
    """
    Calculates the revenue of a divisible item auction based on VCG principles.
    """

    def get_outcome(bids_with_indices):
        """
        Determines the winners and the value each bidder receives in an auction.

        Args:
            bids_with_indices (list of tuples): A list where each tuple is (original_index, bid).

        Returns:
            tuple: A tuple containing:
                - list: A list of winning bidders as (index, bid) tuples.
                - dict: A dictionary mapping each bidder's index to their received value.
        """
        if not bids_with_indices:
            return [], {}

        bids = [b for _, b in bids_with_indices]
        highest_bid = max(bids) if bids else 0
        threshold = 0.5 * highest_bid
        
        winners = [item for item in bids_with_indices if item[1] >= threshold]
        num_winners = len(winners)
        
        outcome_values = {}
        for index, bid in bids_with_indices:
            # If the bidder is a winner, they get an equal share of the item
            if num_winners > 0 and (index, bid) in winners:
                share = 1.0 / num_winners
                outcome_values[index] = bid * share
            else:
                outcome_values[index] = 0
                
        return winners, outcome_values

    def calculate_vcg_revenue(bids):
        """
        Calculates the total auction revenue using the VCG payment rule.

        Args:
            bids (list): A list of bid values.

        Returns:
            float: The total revenue of the auction.
        """
        bids_with_indices = list(enumerate(bids))
        
        # Determine the outcome of the actual auction
        original_winners, original_values = get_outcome(bids_with_indices)
        
        total_revenue = 0
        
        # Calculate payment for each winner based on the harm they cause to others
        for winner_index, winner_bid in original_winners:
            
            # Calculate the total value to others WITH this winner participating
            value_to_others_with_winner = sum(
                val for idx, val in original_values.items() if idx != winner_index
            )

            # Determine the hypothetical outcome WITHOUT this winner
            hypothetical_bids = [
                item for item in bids_with_indices if item[0] != winner_index
            ]
            _, hypothetical_values_for_others = get_outcome(hypothetical_bids)
            
            # The sum of values in the hypothetical scenario is the total value to others
            value_to_others_without_winner = sum(hypothetical_values_for_others.values())
            
            # The winner's payment is the "harm" they caused
            payment = value_to_others_without_winner - value_to_others_with_winner
            total_revenue += payment
            
        return total_revenue

    # Case 1: Bids for calculating x
    bids_x = [100, 20, 5]
    x = calculate_vcg_revenue(bids_x)

    # Case 2: Bids for calculating y
    bids_y = [100, 60, 10]
    y = calculate_vcg_revenue(bids_y)

    print(f"The revenue x for bids {bids_x} is: {int(x)}")
    print(f"The revenue y for bids {bids_y} is: {int(y)}")
    print(f"The final tuple (x, y) is: ({int(x)}, {int(y)})")

solve_divisible_item_auction()