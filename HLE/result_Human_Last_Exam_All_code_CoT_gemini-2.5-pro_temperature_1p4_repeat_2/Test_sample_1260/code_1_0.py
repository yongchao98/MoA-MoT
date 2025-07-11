def calculate_and_explain_returns():
    """
    Demonstrates why buyer's return can be higher than seller's return
    due to differing motivations.
    """
    
    # --- Assumptions ---
    # A stock has a bid price of 99.9 and an ask price of 100.1
    # The initial mid-price is (99.9 + 100.1) / 2 = 100.0
    
    # --- Scenario 1: An "Informed" Buyer ---
    # The buyer believes the price will rise and pays the ask price.
    p_buy_trade = 100.1
    # The trade's information pushes the price up. The new mid-price 15s later is higher.
    p_future_after_buy = 100.2
    
    # The buyer's return is defined as: p'/p - 1
    buyer_return = (p_future_after_buy / p_buy_trade) - 1
    
    print("--- Scenario: Informed Buy ---")
    print(f"Trade Price (p): {p_buy_trade}")
    print(f"Future Mid-Price (p'): {p_future_after_buy}")
    print(f"Buyer's Return Equation: ({p_future_after_buy} / {p_buy_trade}) - 1 = {buyer_return:.6f}")
    print("-" * 20)
    
    # --- Scenario 2: An "Informed" Seller ---
    # The seller believes the price will fall and accepts the bid price.
    p_sell_trade_informed = 99.9
    # The trade's information pushes the price down. The new mid-price 15s later is lower.
    p_future_after_informed_sell = 99.8
    
    # The seller's return is defined as: -(p'/p - 1)
    seller_return_informed = -((p_future_after_informed_sell / p_sell_trade_informed) - 1)

    print("--- Scenario: Informed Sell ---")
    print(f"Trade Price (p): {p_sell_trade_informed}")
    print(f"Future Mid-Price (p'): {p_future_after_informed_sell}")
    print(f"Seller's Return Equation: - (({p_future_after_informed_sell} / {p_sell_trade_informed}) - 1) = {seller_return_informed:.6f}")
    print("-" * 20)

    # --- Scenario 3: An "Uninformed" Seller ---
    # A seller needs cash (liquidity) and accepts the bid price.
    p_sell_trade_uninformed = 99.9
    # This trade carries no new information, so the price impact is temporary and reverts.
    p_future_after_uninformed_sell = 100.0
    
    seller_return_uninformed = -((p_future_after_uninformed_sell / p_sell_trade_uninformed) - 1)

    print("--- Scenario: Uninformed Sell (e.g., for liquidity) ---")
    print(f"Trade Price (p): {p_sell_trade_uninformed}")
    print(f"Future Mid-Price (p'): {p_future_after_uninformed_sell}")
    print(f"Seller's Return Equation: - (({p_future_after_uninformed_sell} / {p_sell_trade_uninformed}) - 1) = {seller_return_uninformed:.6f}")
    print("-" * 20)

    # --- Conclusion: Comparing Average Returns ---
    # Assume buy orders are 100% informed.
    average_buyer_return = buyer_return
    
    # Assume sell orders are 50% informed and 50% uninformed.
    average_seller_return = (seller_return_informed + seller_return_uninformed) / 2
    
    print("\n--- Final Comparison ---")
    print(f"Average Buyer's Return: {average_buyer_return:.6f}")
    print(f"Average Seller's Return: {average_seller_return:.6f}")
    
    is_buyer_return_higher = average_buyer_return > average_seller_return
    print(f"\nIs the Average Buyer's Return higher than the Average Seller's Return? {is_buyer_return_higher}")
    print("\nThis demonstrates that the difference in motivation (informed vs. uninformed/liquidity) can cause the observed asymmetry.")

calculate_and_explain_returns()
<<<L>>>