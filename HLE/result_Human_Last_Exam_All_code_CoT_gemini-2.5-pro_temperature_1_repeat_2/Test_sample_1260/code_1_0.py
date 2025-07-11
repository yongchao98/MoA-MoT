import pandas as pd

def calculate_returns(scenario, mid_price, spread, buy_impact, sell_impact):
    """
    Calculates and prints the 15s returns for a buyer and a seller.

    Args:
        scenario (str): A description of the simulation scenario.
        mid_price (float): The initial mid-price of the stock.
        spread (float): The bid-ask spread.
        buy_impact (float): The price increase after a buy trade.
        sell_impact (float): The price decrease after a sell trade.
    """
    # Initial prices based on spread
    p_ask = mid_price + spread / 2
    p_bid = mid_price - spread / 2

    # --- Buy-Side Order Scenario ---
    # Trade happens at the ask price
    p_trade_buy = p_ask
    # Price moves up due to buy impact
    p_prime_buy = mid_price + buy_impact
    # Calculate buyer's return
    return_buy = (p_prime_buy / p_trade_buy) - 1

    # --- Sell-Side Order Scenario ---
    # Trade happens at the bid price
    p_trade_sell = p_bid
    # Price moves down due to sell impact
    p_prime_sell = mid_price - sell_impact
    # Calculate seller's return
    return_sell = -((p_prime_sell / p_trade_sell) - 1)

    print(f"--- {scenario} ---")
    print(f"Initial Mid-Price: {mid_price:.4f}, Spread: {spread:.4f}")
    print(f"Initial Ask: {p_ask:.4f}, Initial Bid: {p_bid:.4f}\n")

    print("Buyer-Initiated Trade:")
    print(f"  Trade Price (p_buy) = {p_trade_buy:.4f}")
    print(f"  Price 15s later (p'_buy) = {p_prime_buy:.4f}")
    print(f"  Buyer's Return = ({p_prime_buy:.4f} / {p_trade_buy:.4f}) - 1 = {return_buy:.6f} ({return_buy:.4%})")
    print("")

    print("Seller-Initiated Trade:")
    print(f"  Trade Price (p_sell) = {p_trade_sell:.4f}")
    print(f"  Price 15s later (p'_sell) = {p_prime_sell:.4f}")
    print(f"  Seller's Return = - (({p_prime_sell:.4f} / {p_trade_sell:.4f}) - 1) = {return_sell:.6f} ({return_sell:.4%})")
    print("")

    # Compare the returns
    comparison = "higher" if return_buy > return_sell else "lower"
    print(f"Conclusion: Buyer's Return is {comparison} than Seller's Return.\n")
    return return_buy, return_sell

# --- Scenario 1: Symmetric Price Impact ---
# Assume the upward impact of a buy is the same as the downward impact of a sell.
# This tests if the return definition alone causes the phenomenon.
return_buy_sym, return_sell_sym = calculate_returns(
    scenario="Scenario 1: Symmetric Price Impact",
    mid_price=100.0,
    spread=0.02,
    buy_impact=0.03,  # e.g., 3 cents
    sell_impact=0.03   # e.g., 3 cents
)

# --- Scenario 2: Asymmetric Price Impact ---
# Assume the upward impact of a buy is larger than the downward impact of a sell.
# This reflects the theory that buyers' motivations are more information-driven.
return_buy_asym, return_sell_asym = calculate_returns(
    scenario="Scenario 2: Asymmetric Price Impact (Different Motivations)",
    mid_price=100.0,
    spread=0.02,
    buy_impact=0.05,  # e.g., 5 cents
    sell_impact=0.02   # e.g., 2 cents
)

print("="*50)
print("Final Analysis:")
print("Scenario 1 shows that if price impact were symmetric, the seller's return would be higher.")
print("This contradicts the real-world observation.")
print("Scenario 2 shows that when the price impact of buys is larger than for sells, the buyer's return becomes higher.")
print("This matches the observation and shows the key driver is the asymmetric impact of trades, which is best explained by the different motivations of buyers and sellers.")
print("="*50)
