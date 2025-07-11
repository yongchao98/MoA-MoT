# User's information and 2024 IRS limits
user_age = 49
user_turns_50_in_2024 = True

# User's current contributions
user_401k_contribution = 23000
user_ira_contribution = 5000

# 2024 Standard Limits
k401_limit_standard = 23000
ira_limit_standard = 7000

# 2024 Catch-Up Limits (for age 50+)
k401_catch_up = 7500
ira_catch_up = 1000

# --- Calculations ---

# 1. Calculate total 401(k) contribution limit
# Since you turn 50 in 2024, you are eligible for the catch-up contribution.
total_401k_limit = k401_limit_standard + k401_catch_up
remaining_401k = total_401k_limit - user_401k_contribution

print("--- 401(k) Calculation ---")
print(f"Your total 401(k) contribution limit for 2024 (including a ${k401_catch_up} catch-up) is: ${k401_limit_standard} + ${k401_catch_up} = ${total_401k_limit}")
print(f"You have already contributed ${user_401k_contribution}.")
print(f"Remaining 401(k) contribution allowed: ${total_401k_limit} - ${user_401k_contribution} = ${remaining_401k}\n")

# 2. Calculate total IRA contribution limit
# You are also eligible for the IRA catch-up contribution.
total_ira_limit = ira_limit_standard + ira_catch_up
remaining_ira = total_ira_limit - user_ira_contribution

print("--- IRA Calculation ---")
print(f"Your total IRA contribution limit for 2024 (including a ${ira_catch_up} catch-up) is: ${ira_limit_standard} + ${ira_catch_up} = ${total_ira_limit}")
print(f"You have already contributed ${user_ira_contribution}.")
print(f"Remaining IRA contribution allowed: ${total_ira_limit} - ${user_ira_contribution} = ${remaining_ira}\n")

# 3. Calculate total remaining contribution room
total_remaining = remaining_401k + remaining_ira

print("--- Total Remaining Contribution ---")
print("The total additional amount you can legally contribute to your retirement accounts in 2024 is:")
print(f"${remaining_401k} (401k) + ${remaining_ira} (IRA) = ${total_remaining}")

# A note on other accounts
print("\nNote: You have already contributed the maximum amounts to your HSA ($4,150 for self-only coverage) and Health FSA ($3,200) for 2024. These are not included in the retirement calculation above.")

print(f"\n<<<{total_remaining}>>>")