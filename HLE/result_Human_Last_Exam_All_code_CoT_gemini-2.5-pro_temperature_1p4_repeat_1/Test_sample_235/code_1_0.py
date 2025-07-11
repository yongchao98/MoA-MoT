# User's information
user_age_at_year_end = 50 # Turns 50 during 2024, so eligible for catch-up
user_contribution_401k = 23000
user_contribution_ira = 5000

# IRS limits for 2024
limit_401k_base = 23000
limit_401k_catch_up = 7500
limit_ira_base = 7000
limit_ira_catch_up = 1000

# --- 401(k) Calculation ---
print("--- 401(k) Contribution Room ---")
# Calculate total allowed 401k contribution with catch-up
total_allowed_401k = limit_401k_base + limit_401k_catch_up
print(f"Your maximum 401(k) contribution for 2024 (including age 50+ catch-up) is: ${limit_401k_base} + ${limit_401k_catch_up} = ${total_allowed_401k}")

# Calculate remaining 401k contribution room
remaining_401k = total_allowed_401k - user_contribution_401k
print(f"You can contribute an additional: ${total_allowed_401k} (Your Limit) - ${user_contribution_401k} (Your Contribution) = ${remaining_401k}\n")


# --- IRA Calculation ---
print("--- IRA Contribution Room ---")
# Calculate total allowed IRA contribution with catch-up
total_allowed_ira = limit_ira_base + limit_ira_catch_up
print(f"Your maximum IRA contribution for 2024 (including age 50+ catch-up) is: ${limit_ira_base} + ${limit_ira_catch_up} = ${total_allowed_ira}")

# Calculate remaining IRA contribution room
remaining_ira = total_allowed_ira - user_contribution_ira
print(f"You can contribute an additional: ${total_allowed_ira} (Your Limit) - ${user_contribution_ira} (Your Contribution) = ${remaining_ira}\n")


# --- Total Retirement Calculation ---
print("--- Total Remaining Retirement Contributions ---")
# Calculate total remaining contribution room
total_remaining = remaining_401k + remaining_ira
print(f"Total remaining for 2024: ${remaining_401k} (from 401k) + ${remaining_ira} (from IRA) = ${total_remaining}")
print("\nNote: You have already maxed out your HSA ($4,150) and FSA ($3,200) for 2024.")

# Final answer in specified format
# print(f"\n<<<{total_remaining}>>>")
<<<10500>>>