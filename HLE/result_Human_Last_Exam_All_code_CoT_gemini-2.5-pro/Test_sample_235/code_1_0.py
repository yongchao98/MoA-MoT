# User's current status
user_age_turns_in_2024 = 50
user_401k_contribution = 23000
user_hsa_contribution = 4150
user_ira_contribution = 5000

# 2024 IRS Contribution Limits
limit_401k_standard = 23000
limit_401k_catchup = 7500  # For age 50 and over
limit_ira_standard = 7000
limit_ira_catchup = 1000   # For age 50 and over
limit_hsa_self_only = 4150

# --- Calculations ---

# 1. Calculate total 401(k) personal limit for the user
# Since the user turns 50 in 2024, they are eligible for the catch-up.
total_401k_limit = limit_401k_standard + limit_401k_catchup
remaining_401k = total_401k_limit - user_401k_contribution

# 2. Calculate total IRA limit for the user
# Since the user turns 50 in 2024, they are eligible for the catch-up.
total_ira_limit = limit_ira_standard + limit_ira_catchup
remaining_ira = total_ira_limit - user_ira_contribution

# 3. Calculate remaining HSA contribution
remaining_hsa = limit_hsa_self_only - user_hsa_contribution

# 4. Calculate total remaining contributions
total_remaining = remaining_401k + remaining_ira + remaining_hsa

# --- Output the results ---
print("Calculating your remaining retirement contributions for 2024...\n")

print(f"401(k) Calculation:")
print(f"  Your total 401(k) limit (including age 50+ catch-up): ${limit_401k_standard:,.0f} (Standard) + ${limit_401k_catchup:,.0f} (Catch-up) = ${total_401k_limit:,.0f}")
print(f"  You have contributed: ${user_401k_contribution:,.0f}")
print(f"  Remaining 401(k) contribution room: ${remaining_401k:,.0f}\n")

print(f"IRA Calculation:")
print(f"  Your total IRA limit (including age 50+ catch-up): ${limit_ira_standard:,.0f} (Standard) + ${limit_ira_catchup:,.0f} (Catch-up) = ${total_ira_limit:,.0f}")
print(f"  You have contributed: ${user_ira_contribution:,.0f}")
print(f"  Remaining IRA contribution room: ${remaining_ira:,.0f}\n")

print(f"HSA Calculation:")
print(f"  Your HSA limit: ${limit_hsa_self_only:,.0f}")
print(f"  You have contributed: ${user_hsa_contribution:,.0f}")
print(f"  Remaining HSA contribution room: ${remaining_hsa:,.0f}\n")

print("---")
print("Total additional amount you can contribute in 2024:")
print(f"${remaining_401k:,.0f} (401k) + ${remaining_ira:,.0f} (IRA) + ${remaining_hsa:,.0f} (HSA) = ${total_remaining:,.0f}")
<<<10500>>>