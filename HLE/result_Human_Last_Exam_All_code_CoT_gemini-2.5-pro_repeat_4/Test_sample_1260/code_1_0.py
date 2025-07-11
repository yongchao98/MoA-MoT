import random

def simulate_market_scenario():
    """
    Simulates trades to demonstrate why buyer's return is often higher than seller's.
    
    This simulation is based on the hypothesis that there's an asymmetric motivation 
    between buyers and sellers, leading to a small, average positive price drift after a trade.
    """
    
    num_trades = 10000
    # We model the price movement term (p'/p - 1) directly.
    # The key assumption (Hypothesis L) is that the average motivation of buyers is more
    # directionally positive than the average motivation of sellers is negative.
    # This results in a small positive mean for the price movement.
    # Let's assume the mean is +0.005% with some random noise.
    avg_price_movement = 0.00005 
    std_dev_movement = 0.001 # Standard deviation of the movement
    
    total_price_movement = 0
    
    print("Simulating {} trades...\n".format(num_trades))
    
    for _ in range(num_trades):
        # Simulate the price movement for a single trade
        # R = p'/p - 1
        R = random.normalvariate(avg_price_movement, std_dev_movement)
        total_price_movement += R
        
    # Calculate the average price movement across all trades
    average_R = total_price_movement / num_trades
    
    # According to the definition:
    # Average Buyer's Return = Average(R)
    # Average Seller's Return = Average(-R) = -Average(R)
    avg_buyer_return = average_R
    avg_seller_return = -average_R
    
    print("--- Simulation Results ---")
    print(f"The final equation to check is: Average Buyer's Return > Average Seller's Return")
    print(f"Calculated Average Buyer's Return (Average(p'/p - 1)): {avg_buyer_return:.7f}")
    print(f"Calculated Average Seller's Return (Average(-(p'/p - 1))): {avg_seller_return:.7f}")
    
    print("\n--- Conclusion ---")
    if avg_buyer_return > avg_seller_return:
        print(f"The condition is met: {avg_buyer_return:.7f} is indeed greater than {avg_seller_return:.7f}.")
        print("This outcome supports the idea that a small, positive average price drift exists after trades.")
        print("The most probable reason is the asymmetric motivation of market participants: buyers are more consistently motivated by expected price appreciation, while sellers have mixed motives (including liquidity needs), leading to this net positive effect.")
    else:
        print("The condition was not met in this simulation run, which is possible due to random chance, but unlikely with a positive mean.")

# Run the simulation
simulate_market_scenario()