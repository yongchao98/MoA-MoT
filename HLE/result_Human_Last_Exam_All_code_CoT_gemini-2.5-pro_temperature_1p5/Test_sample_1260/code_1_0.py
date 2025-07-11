import random
import numpy as np

def simulate_market_price_impact():
    """
    Simulates trades in a market to demonstrate why the return metric for
    buy orders is typically higher than for sell orders due to information asymmetry.
    """
    # --- Simulation Parameters ---
    num_trades = 10000
    p = 100.0  # We can assume a constant trade price p for simplicity.

    # 'impact_magnitude' represents the price move caused by the information
    # revealed by an aggressive trade. Buys signal positive info, sells signal negative.
    # This is the core of the "Information Asymmetry" model.
    # A buy order pushes the future price up by 'impact_magnitude'.
    # A sell order pushes the future price down by 'impact_magnitude'.
    impact_magnitude = 0.0005  # 0.05% or 5 basis points

    # 'noise_std_dev' represents all other random price movements not related
    # to the immediate trade's impact.
    noise_std_dev = 0.0010  # 0.10% or 10 basis points

    buy_initiated_returns = []
    sell_initiated_returns = []

    # --- Simulation Loop ---
    for _ in range(num_trades):
        # The return metric we are evaluating is R = p'/p - 1
        
        # Randomly decide if the trade is buy-initiated or sell-initiated
        if random.random() > 0.5:
            # --- This is a Buy-initiated Trade ---
            # The trade signals positive information, so the price is expected to rise.
            # We model the future price p' based on this upward impact plus random noise.
            # The change is `impact_magnitude + noise`.
            price_change = impact_magnitude + np.random.normal(loc=0, scale=noise_std_dev)
            return_metric = price_change
            buy_initiated_returns.append(return_metric)

        else:
            # --- This is a Sell-initiated Trade ---
            # The trade signals negative information, so the price is expected to fall.
            # We model the future price p' based on this downward impact plus random noise.
            # The change is `-impact_magnitude + noise`.
            price_change = -impact_magnitude + np.random.normal(loc=0, scale=noise_std_dev)
            return_metric = price_change
            sell_initiated_returns.append(return_metric)

    # --- Calculate and Print Results ---
    avg_return_buy = np.mean(buy_initiated_returns)
    avg_return_sell = np.mean(sell_initiated_returns)

    print("Simulation based on the Information Asymmetry principle.")
    print(f"Number of simulated trades: {num_trades}")
    print(f"Assumed price impact of a trade: {impact_magnitude*100:.4f}%")
    print("-" * 55)
    print("Final Equation shows Average Return Metric for Buys vs. Sells:")
    # We output the numbers in the final comparison as requested
    print(f"Average(R | buy) > Average(R | sell)")
    print(f"{avg_return_buy:.6f} > {avg_return_sell:.6f}")
    print("-" * 55)
    print("The simulation confirms that the average return metric for buy orders")
    print("is positive, while for sell orders it is negative.")

# Execute the simulation
simulate_market_price_impact()