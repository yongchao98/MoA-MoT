# 2024 Tax Year Contribution Analysis

# User Information
user_age_turns_in_2024 = 50
user_401k_contribution = 23000
user_ira_contribution = 5000
user_hsa_contribution = 4150
user_fsa_contribution = 3200

# IRS Limits for 2024
# Note: The user turns 50 in 2024, so they are eligible for catch-up contributions for the entire year.
standard_401k_limit = 23000
catch_up_401k = 7500
total_401k_limit = standard_401k_limit + catch_up_401k

standard_ira_limit = 7000
catch_up_ira = 1000
total_ira_limit = standard_ira_limit + catch_up_ira

hsa_limit_self = 4150
fsa_limit = 3200

# --- Calculations ---

# 1. 401k Calculation
remaining_401k = total_401k_limit - user_401k_contribution
print("--- 401k Contribution Room ---")
print(f"Your Total 401k Contribution Limit (Standard + Age 50+ Catch-up): ${total_401k_limit:,.2f}")
print(f"Your Contributions So Far: ${user_401k_contribution:,.2f}")
print(f"Remaining 401k Contribution for 2024: ${remaining_401k:,.2f}\n")

# 2. IRA Calculation (for Backdoor Roth)
remaining_ira = total_ira_limit - user_ira_contribution
print("--- IRA Contribution Room ---")
print(f"Your Total IRA Contribution Limit (Standard + Age 50+ Catch-up): ${total_ira_limit:,.2f}")
print(f"Your Contributions So Far: ${user_ira_contribution:,.2f}")
print(f"Remaining IRA Contribution for 2024: ${remaining_ira:,.2f}\n")

# 3. Other Accounts (for context)
print("--- Other Accounts ---")
print(f"Your HSA contribution of ${user_hsa_contribution:,.2f} has reached the annual limit of ${hsa_limit_self:,.2f}.")
print("Note: An FSA is not a retirement account.\n")

# 4. Final Summary
total_remaining_contribution = remaining_401k + remaining_ira
print("--- Total Remaining Retirement Contributions for 2024 ---")
print("This is the total additional amount you can contribute across your retirement accounts.")
print(f"${remaining_401k:,.2f} (from 401k) + ${remaining_ira:,.2f} (from IRA) = ${total_remaining_contribution:,.2f}")

<<<10500.00>>>