def analyze_returns(trade_price, future_mid_price):
    """
    This function demonstrates the relationship between buyer and seller returns.
    It shows that the buyer's return is higher than the seller's return
    if and only if the future mid-price is greater than the trade price.
    """
    p = trade_price
    p_prime = future_mid_price

    print(f"--- Analyzing Scenario: Trade Price (p) = {p:.2f}, Future Mid-Price (p') = {p_prime:.2f} ---")

    # Define the core component of the return calculation
    base_return = (p_prime / p) - 1

    # Calculate returns for buyer and seller
    buyer_return = base_return
    seller_return = -base_return

    # Print the full equations as requested
    print("Buyer's Return Calculation:")
    print(f"  R_buy = p' / p - 1")
    print(f"  R_buy = {p_prime:.2f} / {p:.2f} - 1 = {buyer_return:.6f}")
    print(f"  R_buy = {buyer_return:.4%}\n")

    print("Seller's Return Calculation:")
    print(f"  R_sell = -(p' / p - 1)")
    print(f"  R_sell = -({p_prime:.2f} / {p:.2f} - 1) = {seller_return:.6f}")
    print(f"  R_sell = {seller_return:.4%}\n")

    # Check the condition from the problem
    is_higher = buyer_return > seller_return
    print(f"Condition Check: Is Buyer's Return > Seller's Return?")
    print(f"Result: {buyer_return:.4%} > {seller_return:.4%}, which is {is_higher}.")
    print("-" * 70)


print("The problem states that for a trade, the 15s return is usually higher for buyers than for sellers.")
print("This script will demonstrate the mathematical condition required for this to be true.\n")

# Case 1: The scenario described by Alice, where the future price is higher.
# This aligns with the reasoning that price tends to increase after a trade on average.
print("Case 1: Future price is HIGHER than trade price (p' > p)")
analyze_returns(trade_price=100.00, future_mid_price=100.04)

# Case 2: The opposite scenario, where the future price is lower.
print("Case 2: Future price is LOWER than trade price (p' < p)")
analyze_returns(trade_price=100.00, future_mid_price=99.96)

print("Conclusion from the code:")
print("The buyer's return is higher than the seller's return only in Case 1, where p' > p.")
print("Thus, Alice's observation implies that the market exhibits a tendency for p' > p, and the reasoning provided above explains why this occurs.")
