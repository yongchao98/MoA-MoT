import pandas as pd

def analyze_market_scenario():
    """
    This function simulates a few trades to demonstrate why the average 15s return
    for buyers can be higher than for sellers, based on the principle of
    information asymmetry.
    """
    # Let's create a few hypothetical trades.
    # A trade is initiated by a trader who has information (or a strong belief)
    # about the future price direction.
    # Trade Type: The initiator of the trade (Buyer or Seller).
    # p: The price at which the trade occurred.
    # p_prime: The mid-price 15 seconds after the trade.
    
    trades_data = [
        {'trade_id': 1, 'type': 'Buy-initiated', 'p': 100.10, 'p_prime': 100.15}, # Informed buyer pushes price up
        {'trade_id': 2, 'type': 'Sell-initiated', 'p': 100.05, 'p_prime': 100.00}, # Informed seller pushes price down
        {'trade_id': 3, 'type': 'Buy-initiated', 'p': 100.20, 'p_prime': 100.28}, # Another informed buyer
        {'trade_id': 4, 'type': 'Buy-initiated', 'p': 100.30, 'p_prime': 100.36}, # Markets often have more buy pressure
    ]
    
    df = pd.DataFrame(trades_data)
    
    print("--- Hypothetical Trade Data ---")
    print(df)
    print("\n--- Calculating 15s Returns ---")

    # The problem defines returns for the buyer and seller in *any* given trade.
    # Let's calculate R = p'/p - 1 for each trade first.
    # Buyer's Return = R
    # Seller's Return = -R
    
    buyer_returns = []
    seller_returns = []

    for index, row in df.iterrows():
        p = row['p']
        p_prime = row['p_prime']
        
        # Calculate R = p'/p - 1
        R = (p_prime / p) - 1
        
        # This R is the buyer's return
        buyer_return = R
        buyer_returns.append(buyer_return)
        
        # The seller's return is -R
        seller_return = -R
        seller_returns.append(seller_return)
        
        print(f"\nFor Trade {row['trade_id']} (p={p:.4f}, p'={p_prime:.4f}):")
        print(f"  The term (p'/p - 1) is: ({p_prime:.4f} / {p:.4f}) - 1 = {R:.6f}")
        print(f"  Buyer's 15s Return = {buyer_return:.6f}")
        print(f"  Seller's 15s Return = {seller_return:.6f}")

    df['buyer_return'] = buyer_returns
    df['seller_return'] = seller_returns

    # Calculate the average returns across all trades
    avg_buyer_return = df['buyer_return'].mean()
    avg_seller_return = df['seller_return'].mean()

    print("\n--- Final Analysis ---")
    print(f"Average Buyer's Return: {avg_buyer_return:.6f}")
    print(f"Average Seller's Return: {avg_seller_return:.6f}")
    
    # Check the condition from the problem
    is_higher = avg_buyer_return > avg_seller_return
    
    print(f"\nIs the Average Buyer's Return > Average Seller's Return? {is_higher}")
    
    # This is equivalent to checking if Average(R) > 0
    avg_R = (df['buyer_return']).mean()
    print(f"This is because the average of (p'/p - 1) is {avg_R:.6f}, which is greater than 0.")
    print("\nThis positive average is explained by information asymmetry, where informed trades move the price and a general market updrift means buy-initiated trades and their effects are more common or pronounced.")

# Execute the analysis
analyze_market_scenario()