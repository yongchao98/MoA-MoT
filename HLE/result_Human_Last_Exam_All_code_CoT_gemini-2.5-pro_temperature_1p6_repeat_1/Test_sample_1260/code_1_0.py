# Plan:
# 1. Define initial parameters for a hypothetical trade.
# 2. Model the price impact, making the impact of a buy trade stronger
#    than that of a sell trade, reflecting the "different motivations" theory.
# 3. Calculate the buyer's and seller's 15s return based on the problem's definition.
# 4. Print the step-by-step calculation to show that the buyer's return is indeed higher.

# Step 1: Define initial parameters
p = 100.0  # Initial trade price

print(f"Let's assume a trade happens at a price (p) of: ${p:.2f}\n")

# Step 2: Model the price impact
# The core hypothesis: Buyers are more often motivated by information, leading to a stronger price impact.
# Sellers are often motivated by liquidity needs, leading to a weaker price impact.
# Let's say a buy trade leads to a 0.10% price increase after 15s.
# And a sell trade leads to a 0.05% price decrease after 15s (smaller magnitude).

# Price 15s after a BUY trade
p_prime_after_buy = 100.10
# Price 15s after a SELL trade
p_prime_after_sell = 99.95

print("--- Calculating Buyer's 15s Return ---")
print(f"Due to strong positive information, the price 15s after a BUY (p') is: ${p_prime_after_buy:.2f}")
# Step 3: Calculate buyer's return
buyer_return = (p_prime_after_buy / p) - 1
print(f"Buyer's Return = (p' / p) - 1")
print(f"Buyer's Return = ({p_prime_after_buy:.2f} / {p:.2f}) - 1 = {buyer_return:.4f}\n")


print("--- Calculating Seller's 15s Return ---")
print(f"Due to weaker negative information, the price 15s after a SELL (p') is: ${p_prime_after_sell:.2f}")
# Step 4: Calculate seller's return
seller_return = -((p_prime_after_sell / p) - 1)
print(f"Seller's Return = -((p' / p) - 1)")
print(f"Seller's Return = -(({p_prime_after_sell:.2f} / {p:.2f}) - 1) = -({p_prime_after_sell/p - 1:.4f}) = {seller_return:.4f}\n")

# Step 5: Compare the returns
print("--- Comparison ---")
print(f"Buyer's Return ({buyer_return:.4f}) vs. Seller's Return ({seller_return:.4f})")

is_higher = buyer_return > seller_return
print(f"Is the buyer's return higher than the seller's return? {is_higher}")
print(f"The calculation shows: {buyer_return:.4f} > {seller_return:.4f}")
