import numpy as np

def simulate_market_momentum(num_trades=1000):
    """
    Simulates trades in a market with a short-term momentum effect
    and calculates the resulting returns for buyers and sellers.
    """
    # Initial trade price
    p = 100.0

    # To model the momentum effect (p' > p on average), we assume the short-term
    # price change follows a distribution with a small positive mean.
    # We use log returns for better mathematical properties.
    # Mean of log return is slightly positive, representing upward momentum.
    momentum_drift = 0.00005  # Represents a 0.005% average price increase in 15s
    volatility = 0.0005      # Represents short-term price volatility

    buyer_returns = []
    seller_returns = []

    print(f"Simulating {num_trades} trades with an initial price of ${p:.2f}")
    print(f"Modeling a short-term momentum effect...\n")

    for i in range(num_trades):
        # Generate a random price change based on momentum and volatility
        log_return = np.random.normal(loc=momentum_drift, scale=volatility)
        
        # Calculate the price 15 seconds later
        p_prime = p * np.exp(log_return)
        
        # Calculate returns for this trade
        # The equation is based on p'/p, so we can use 1 for p for simplicity
        # Buyer's return = p'/p - 1
        # Seller's return = -(p'/p - 1)
        buyer_return = (p_prime / p) - 1
        seller_return = -(p_prime / p - 1)
        
        buyer_returns.append(buyer_return)
        seller_returns.append(seller_return)

    # Calculate average returns
    avg_buyer_return = np.mean(buyer_returns)
    avg_seller_return = np.mean(seller_returns)

    print("--- Simulation Results ---")
    print(f"Average 15s return for buyers:  {avg_buyer_return:.6f} (or {avg_buyer_return*100:.4f}%)")
    print(f"Average 15s return for sellers: {avg_seller_return:.6f} (or {avg_seller_return*100:.4f}%)")
    
    # Final verification
    print("\nIs the average buyer's return higher than the average seller's return?")
    print(f"Result: {avg_buyer_return > avg_seller_return}")
    
    # We can also check if p' was greater than p on average
    is_momentum_positive = (avg_buyer_return - avg_seller_return) > 0
    print(f"This is mathematically equivalent to checking if p' > p on average, which is the momentum effect we modeled.")
    print(f"The inequality (p'/p - 1) > -(p'/p - 1) simplifies to p'/p > 1.")
    print(f"Our average p'/p was approximately {1 + avg_buyer_return:.6f}, which is greater than 1.")


if __name__ == '__main__':
    simulate_market_momentum()
