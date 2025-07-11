# User's current contributions
current_401k_contribution = 23000
current_ira_contribution = 5000
user_age_in_2024 = 50

# 2024 IRS Contribution Limits
limit_401k_base = 23000
limit_401k_catch_up = 7500
limit_ira_base = 7000
limit_ira_catch_up = 1000

# Since the user turns 50 in 2024, they are eligible for catch-up contributions.
total_401k_limit = limit_401k_base + limit_401k_catch_up
total_ira_limit = limit_ira_base + limit_ira_catch_up

# --- 401(k) Calculation ---
print("--- 401(k) Contribution ---")
print(f"Your total 401(k) contribution limit for 2024 (including age 50+ catch-up) is: ${limit_401k_base} + ${limit_401k_catch_up} = ${total_401k_limit}")
print(f"You have already contributed: ${current_401k_contribution}")
remaining_401k = total_401k_limit - current_401k_contribution
print(f"Remaining 401(k) contribution allowed: ${total_401k_limit} - ${current_401k_contribution} = ${remaining_401k}\n")


# --- IRA Contribution Calculation ---
print("--- IRA Contribution ---")
print(f"Your total IRA contribution limit for 2024 (including age 50+ catch-up) is: ${limit_ira_base} + ${limit_ira_catch_up} = ${total_ira_limit}")
print(f"You have already contributed: ${current_ira_contribution}")
remaining_ira = total_ira_limit - current_ira_contribution
print(f"Remaining IRA contribution allowed: ${total_ira_limit} - ${current_ira_contribution} = ${remaining_ira}\n")

# --- Total Calculation ---
print("--- Total Remaining Contribution for 2024 ---")
total_remaining = remaining_401k + remaining_ira
print(f"Your HSA and FSA are already maxed out for 2024.")
print(f"Total additional amount you can contribute to retirement accounts is:")
print(f"${remaining_401k} (from 401k) + ${remaining_ira} (from IRA) = ${total_remaining}")
