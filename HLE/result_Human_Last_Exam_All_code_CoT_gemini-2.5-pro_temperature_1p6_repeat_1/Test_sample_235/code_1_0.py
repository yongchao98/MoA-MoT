# 2024 IRS Contribution Limits
limit_401k_standard = 23000
limit_ira_standard = 7000

# Age 50+ Catch-up Contributions (eligible since user turns 50 in 2024)
catch_up_401k = 7500
catch_up_ira = 1000

# User's Current Contributions
contrib_401k = 23000
contrib_ira = 5000

# Other tax-advantaged accounts (already maxed out for 2024)
# HSA Limit (Self-Only): $4,150. User contributed $4,150.
# FSA Limit: $3,200. User contributed $3,200.

# 1. Calculate total allowable 401(k) contribution for the user
total_limit_401k = limit_401k_standard + catch_up_401k

# 2. Calculate remaining 401(k) contribution room
remaining_401k = total_limit_401k - contrib_401k

# 3. Calculate total allowable IRA contribution for the user
total_limit_ira = limit_ira_standard + catch_up_ira

# 4. Calculate remaining IRA contribution room
remaining_ira = total_limit_ira - contrib_ira

# 5. Calculate total remaining retirement contributions
total_remaining = remaining_401k + remaining_ira

# Print the breakdown and final answer
print("Calculating your remaining retirement contribution room for 2024...")
print("-" * 60)
print(f"Your total 401(k) contribution limit (including ${catch_up_401k:,} catch-up): ${total_limit_401k:,}")
print(f"Your total IRA contribution limit (including ${catch_up_ira:,} catch-up): ${total_limit_ira:,}")
print("\nYour Health Savings Account (HSA) and Flexible Spending Account (FSA) are fully funded for 2024.")
print("-" * 60)
print("\nHere is the calculation for the additional amount you can contribute:\n")

# Final Equation Output
print(f"Remaining 401(k) Contribution Room: (${total_limit_401k:,} limit - ${contrib_401k:,} contributed) = ${remaining_401k:,}")
print(f"Remaining IRA Contribution Room:   (${total_limit_ira:,} limit - ${contrib_ira:,} contributed) = ${remaining_ira:,}")
print(f"\nTotal Remaining Contribution = ${remaining_401k:,} (401k) + ${remaining_ira:,} (IRA) = ${total_remaining:,}")
