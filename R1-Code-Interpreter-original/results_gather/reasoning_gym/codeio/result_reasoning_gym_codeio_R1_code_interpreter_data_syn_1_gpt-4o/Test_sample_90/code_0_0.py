import math
import collections

def calculate_scores(market_name, house_bid, player_bids, outcome):
    # Initialize the market with the house bid
    bids = [('House', house_bid)]
    results = collections.defaultdict(float)
    
    # Add player bids
    for player, bid in player_bids:
        bids.append((player, bid))
    
    # Resolve the market
    if not outcome:
        # If the outcome is False, use the complement of the bids
        bids = [(player, 1 - bid) for player, bid in bids]
    
    # Calculate scores
    _, last_bid = bids[0]
    for player, bid in bids[1:]:
        log_prob = math.log(bid, 2)
        log_prob_last = math.log(last_bid, 2)
        score = 100 * (log_prob - log_prob_last)
        results[player] += score
        last_bid = bid
    
    return dict(results)

# Given input
market_name = 'Market_28'
house_bid = 0.5969980871075791
player_bids = [('Player_1', 0.9095422098787996), ('Player_2', 0.07098421998126647), ('Player_3', 0.36655074924986386)]
outcome = True

# Calculate and print the scores
scores = calculate_scores(market_name, house_bid, player_bids, outcome)
print(scores)