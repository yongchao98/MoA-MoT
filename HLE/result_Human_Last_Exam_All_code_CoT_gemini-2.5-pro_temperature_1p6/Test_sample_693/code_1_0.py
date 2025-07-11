def calculate_auction_revenue(bids):
    """
    Calculates the revenue of a divisible item auction based on the given rules.

    The auction is truthful, and the item is equally divided among bidders whose
    bid is at least 1/2 of the highest bid. The revenue is determined by the
    critical value set by the highest losing bid.

    Revenue = 2 * highest_losing_bid
    """
    if not bids or len(bids) < 2:
        return 0, 0, 0, [] # revenue, b_max, threshold, highest_loser_bid

    b_max = max(bids)
    threshold = 0.5 * b_max

    winners = []
    losers = []
    for bid in bids:
        if bid >= threshold:
            winners.append(bid)
        else:
            losers.append(bid)

    if not losers:
        return 0, b_max, threshold, None

    highest_losing_bid = max(losers)
    revenue = 2 * highest_losing_bid
    return revenue, b_max, threshold, highest_losing_bid

# Bids for the two scenarios
bids_x = [100, 20, 5]
bids_y = [100, 60, 10]

# Calculate revenue for x
x, b_max_x, threshold_x, highest_loser_x = calculate_auction_revenue(bids_x)

# Calculate revenue for y
y, b_max_y, threshold_y, highest_loser_y = calculate_auction_revenue(bids_y)

# Print the results
print("Calculating the revenue for each scenario:")

# Scenario x
print("\nScenario for x with bids " + str(bids_x) + ":")
print(f"1. The highest bid is {b_max_x}.")
print(f"2. The winning threshold is 0.5 * {b_max_x} = {threshold_x}.")
print(f"3. The losers are bidders with bids less than {threshold_x}. The highest losing bid is {highest_loser_x}.")
print(f"4. The revenue x is 2 * {highest_loser_x} = {x}.")


# Scenario y
print("\nScenario for y with bids " + str(bids_y) + ":")
print(f"1. The highest bid is {b_max_y}.")
print(f"2. The winning threshold is 0.5 * {b_max_y} = {threshold_y}.")
print(f"3. The losers are bidders with bids less than {threshold_y}. The highest losing bid is {highest_loser_y}.")
print(f"4. The revenue y is 2 * {highest_loser_y} = {y}.")


print(f"\nFinal result (x, y) = ({x}, {y})")
