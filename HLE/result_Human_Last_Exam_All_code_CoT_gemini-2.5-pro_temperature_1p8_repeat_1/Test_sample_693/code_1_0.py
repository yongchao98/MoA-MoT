def calculate_auction_revenue(bids, case_variable):
    """
    Calculates the revenue from a divisible item auction based on the specified rules.

    Args:
        bids (list): A list of numerical bids.
        case_variable (str): The variable name for the case ('x' or 'y') for printing.
    """
    print(f"Calculating revenue for case '{case_variable}' with bids: {bids}")

    if not bids:
        print("No bids provided. Revenue is 0.")
        return 0

    highest_bid = max(bids)
    print(f"The highest bid is: {highest_bid}")

    winning_threshold = 0.5 * highest_bid
    print(f"The winning threshold is 0.5 * {highest_bid} = {winning_threshold}")

    winners = [b for b in bids if b >= winning_threshold]
    losers = sorted([b for b in bids if b < winning_threshold], reverse=True)

    num_winners = len(winners)
    print(f"The winning bids are: {winners} ({num_winners} winner(s))")

    if not losers:
        # This case happens if all bidders win
        price = 0
        print("All bidders won, so the price is 0.")
    else:
        price = losers[0]
        print(f"The losing bids are: {sorted(losers, reverse=True)}")
        print(f"The price is the highest losing bid: {price}")

    revenue = num_winners * price
    # The final equation is requested in the output format
    print(f"The revenue calculation is: {num_winners} (winners) * {price} (price)")
    print(f"Final revenue for '{case_variable}' is: {revenue}\n")
    return revenue

def main():
    """
    Solves the auction problem for the two given cases and prints the result.
    """
    # Case 1 for x
    bids_x = [100, 20, 5]
    x = calculate_auction_revenue(bids_x, 'x')

    # Case 2 for y
    bids_y = [100, 60, 10]
    y = calculate_auction_revenue(bids_y, 'y')

    # Final result
    print("The final result (x, y) is:")
    print(f"({x}, {y})")

if __name__ == "__main__":
    main()
<<< (20, 20) >>>