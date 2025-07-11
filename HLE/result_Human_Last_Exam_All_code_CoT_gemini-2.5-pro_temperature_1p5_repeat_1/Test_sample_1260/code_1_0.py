# --- Plan ---
# 1. Define a hypothetical trade price (p) and a mid-price 15 seconds later (p').
# 2. To align with the scenario, p' will be slightly higher than p.
# 3. Calculate the buyer's 15s return using the formula: p'/p - 1.
# 4. Calculate the seller's 15s return using the formula: -(p'/p - 1).
# 5. Print the step-by-step calculations for both returns to show how the numbers are derived.
# 6. Conclude by showing that the buyer's return is greater than the seller's in this case.

# --- Code Execution ---

# Define the trade prices
p = 100.00 # Trade price
p_prime = 100.05 # Mid-price 15 seconds later

# Calculate the returns
buyer_return = (p_prime / p) - 1
seller_return = -((p_prime / p) - 1)

print("This code demonstrates the calculation for a single trade.")
print("The key observation is that the buyer's 15s return is usually higher than the seller's.")
print("This happens when the price 15 seconds later (p') is higher than the trade price (p).")
print("-" * 30)

print(f"Let's assume a trade price (p) = {p}")
print(f"Let's assume the mid-price 15s later (p') = {p_prime}")
print("-" * 30)

# Display the buyer's return calculation
print("Buyer's 15s Return = p' / p - 1")
print(f"                     = {p_prime} / {p} - 1")
print(f"                     = {p_prime / p:.4f} - 1")
print(f"                     = {buyer_return:.4f}")
print("")

# Display the seller's return calculation
print("Seller's 15s Return = -(p' / p - 1)")
print(f"                     = -({p_prime} / {p} - 1)")
print(f"                     = -({p_prime / p:.4f} - 1)")
print(f"                     = -({buyer_return:.4f})")
print(f"                     = {seller_return:.4f}")
print("-" * 30)

print(f"Result: The buyer's return ({buyer_return:.4f}) is higher than the seller's return ({seller_return:.4f}).")
print("This asymmetry is likely caused by systematic differences in motivation and urgency between buyers and sellers.")
<<<L>>>