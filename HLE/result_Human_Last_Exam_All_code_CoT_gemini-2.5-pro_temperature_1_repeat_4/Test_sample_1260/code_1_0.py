import sys

# Plan:
# 1. Set up an example market state with a bid and an ask price.
# 2. Define the trade prices for a buyer (p_buy) and a seller (p_sell).
# 3. Assume a future mid-price (p_prime). To isolate the core effect, we'll
#    assume the price is stable.
# 4. Calculate the returns for the buyer and seller based on the provided definitions.
# 5. Show the step-by-step algebraic proof that the buyer's return is structurally
#    higher than the seller's, simplifying to the expression (ask - bid)^2 > 0.
# 6. Print the final equation with the example numbers to confirm the logic.

# Step 1: Set up example market prices
bid_price = 10.00
ask_price = 10.02 # In any market, ask > bid

# Step 2: A buyer's trade price 'p' is the ask, a seller's is the bid.
p_buy = ask_price
p_sell = bid_price

# Step 3: Let's assume the mid-price 15 seconds later, p_prime, is the same as the
# initial mid-price. This isolates the effect of the bid-ask spread itself.
p_prime = (ask_price + bid_price) / 2

# Step 4: Calculate the returns based on the definitions
# Buyer's 15s return = p'/p - 1
r_buy = (p_prime / p_buy) - 1
# Seller's 15s return = -(p'/p - 1)
r_sell = -((p_prime / p_sell) - 1)

print("--- Numerical Example ---")
print(f"Buyer's trade price (ask): {p_buy}")
print(f"Seller's trade price (bid): {p_sell}")
print(f"Assumed future mid-price (p_prime): {p_prime:.4f}")
print(f"Buyer's Return = ({p_prime:.4f} / {p_buy}) - 1 = {r_buy:.6f}")
print(f"Seller's Return = -(({p_prime:.4f} / {p_sell}) - 1) = {r_sell:.6f}")
print(f"Result: The buyer's return ({r_buy:.6f}) is indeed higher than the seller's return ({r_sell:.6f}).")

print("\n--- Why This Happens: The Algebraic Proof ---")
print("The observation is a direct result of the return definitions and the bid-ask spread.")
print("The inequality to test is: Buyer's Return > Seller's Return")
print("1. (p_prime / ask_price) - 1 > -((p_prime / bid_price) - 1)")
print("2. (p_prime / ask_price) - 1 > 1 - (p_prime / bid_price)")
print("3. (p_prime / ask_price) + (p_prime / bid_price) > 2")
print("4. p_prime * ((ask_price + bid_price) / (ask_price * bid_price)) > 2")
print("5. Assuming p_prime is the mid-price ((ask_price + bid_price) / 2):")
print("   ((ask_price + bid_price) / 2) * ((ask_price + bid_price) / (ask_price * bid_price)) > 2")
print("6. (ask_price + bid_price)^2 / (2 * ask_price * bid_price) > 2")
print("7. (ask_price + bid_price)^2 > 4 * ask_price * bid_price")
print("8. ask_price^2 + 2*ask_price*bid_price + bid_price^2 > 4 * ask_price * bid_price")
print("9. ask_price^2 - 2*ask_price*bid_price + bid_price^2 > 0")
print("10. (ask_price - bid_price)^2 > 0")

print("\n--- Final Equation Check ---")
print("The final derived equation is (ask_price - bid_price)^2 > 0.")
print("This is always true as long as the bid-ask spread (ask_price - bid_price) is not zero.")
print("Plugging in our example numbers for each part of the final equation:")
lhs_val = (ask_price - bid_price)**2
rhs_val = 0
print(f"Equation: ({ask_price} - {bid_price})^2 > {rhs_val}")
print(f"Result: {lhs_val:.4f} > {rhs_val}")
print(f"This is {lhs_val > rhs_val}, confirming the logic.")
print("\nSince all major exchanges (SZSE, SSE, Nasdaq, HKEX) have a bid-ask spread, this phenomenon is expected in all of them. The reason is the nature of the return definition.")

<<<A>>>