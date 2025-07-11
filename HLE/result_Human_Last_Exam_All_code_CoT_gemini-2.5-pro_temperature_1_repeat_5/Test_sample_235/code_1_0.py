# User's current financial information
age = 49
birthday_month = 6 # June
user_401k_contribution = 23000
user_ira_contribution = 5000

# IRS limits for 2024
# The user turns 50 in 2024, making them eligible for catch-up contributions for the entire year.
standard_401k_limit = 23000
catch_up_401k = 7500
total_401k_limit = standard_401k_limit + catch_up_401k

standard_ira_limit = 7000
catch_up_ira = 1000
total_ira_limit = standard_ira_limit + catch_up_ira

# --- Calculations ---

# 1. Calculate remaining 401(k) contribution room
remaining_401k = total_401k_limit - user_401k_contribution
print(f"Your total 401(k) contribution limit for 2024 (including age 50+ catch-up) is ${total_401k_limit:,.2f}.")
print(f"Remaining 401(k) contribution = ${total_401k_limit:,.0f} (Total Limit) - ${user_401k_contribution:,.0f} (Your Contribution) = ${remaining_401k:,.2f}")
print("-" * 30)

# 2. Calculate remaining IRA contribution room
remaining_ira = total_ira_limit - user_ira_contribution
print(f"Your total IRA contribution limit for 2024 (including age 50+ catch-up) is ${total_ira_limit:,.2f}.")
print(f"Remaining IRA contribution = ${total_ira_limit:,.0f} (Total Limit) - ${user_ira_contribution:,.0f} (Your Contribution) = ${remaining_ira:,.2f}")
print("-" * 30)

# 3. Calculate total remaining retirement contributions
total_remaining_contribution = remaining_401k + remaining_ira
print("Total additional amount you can contribute to retirement accounts in 2024:")
print(f"${remaining_401k:,.0f} (from 401k) + ${remaining_ira:,.0f} (from IRA) = ${total_remaining_contribution:,.2f}")
