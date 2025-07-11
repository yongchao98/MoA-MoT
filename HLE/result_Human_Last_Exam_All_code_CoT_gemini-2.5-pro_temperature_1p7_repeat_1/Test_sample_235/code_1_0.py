# User's financial information for 2024
user_age = 49
user_turns_50_in_2024 = True # Birthday is June 3rd
contribution_401k_current = 23000
contribution_ira_current = 5000

# 2024 IRS Contribution Limits
limit_401k_standard = 23000
limit_401k_catch_up = 7500
limit_ira_standard = 7000
limit_ira_catch_up = 1000

# --- 401(k) Calculation ---
print("--- 401(k) Contribution Analysis ---")
total_401k_limit = limit_401k_standard
if user_turns_50_in_2024:
    total_401k_limit += limit_401k_catch_up

remaining_401k = total_401k_limit - contribution_401k_current
print(f"Your total 401(k) contribution limit for 2024 is ${total_401k_limit:,} (${limit_401k_standard:,} standard + ${limit_401k_catch_up:,} catch-up).")
print(f"You have contributed ${contribution_401k_current:,}, leaving you with ${remaining_401k:,} remaining to contribute.\n")


# --- IRA Calculation ---
print("--- IRA Contribution Analysis ---")
total_ira_limit = limit_ira_standard
if user_turns_50_in_2024:
    total_ira_limit += limit_ira_catch_up

remaining_ira = total_ira_limit - contribution_ira_current
print(f"Your total IRA contribution limit for 2024 is ${total_ira_limit:,} (${limit_ira_standard:,} standard + ${limit_ira_catch_up:,} catch-up).")
print(f"You have contributed ${contribution_ira_current:,}, leaving you with ${remaining_ira:,} remaining to contribute.\n")

# --- Total Remaining Contribution Calculation ---
print("--- Total Remaining Retirement Contribution ---")
total_remaining = remaining_401k + remaining_ira
print(f"The total additional amount you can contribute to your retirement accounts in 2024 is:")
print(f"${remaining_401k:,} (401k) + ${remaining_ira:,} (IRA) = ${total_remaining:,}")

<<<10500>>>