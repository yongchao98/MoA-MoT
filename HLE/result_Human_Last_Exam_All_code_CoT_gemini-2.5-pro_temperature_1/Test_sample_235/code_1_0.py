# User's current situation
user_age = 49
user_turns_50_in_2024 = True
user_401k_contribution = 23000
user_ira_contribution = 5000
user_hsa_contribution = 4150

# 2024 IRS Contribution Limits
limit_401k_standard = 23000
limit_401k_catchup = 7500  # For age 50+
limit_ira_standard = 7000
limit_ira_catchup = 1000   # For age 50+
limit_hsa_self_only = 4150

# --- 1. Calculate 401(k) Room ---
# Since the user turns 50 in 2024, they are eligible for the catch-up contribution.
total_allowed_401k = limit_401k_standard + limit_401k_catchup
remaining_401k = total_allowed_401k - user_401k_contribution

print("--- 401(k) Calculation ---")
print(f"Your total allowed 401(k) contribution for 2024 (including a ${limit_401k_catchup} catch-up) is: ${total_allowed_401k}")
print(f"You have already contributed: ${user_401k_contribution}")
print(f"Equation: ${total_allowed_401k} (Total Allowed) - ${user_401k_contribution} (Contributed) = ${remaining_401k}")
print(f"Remaining 401(k) contribution room: ${remaining_401k}\n")

# --- 2. Calculate IRA Room ---
# Since the user turns 50 in 2024, they are also eligible for the IRA catch-up.
total_allowed_ira = limit_ira_standard + limit_ira_catchup
remaining_ira = total_allowed_ira - user_ira_contribution

print("--- IRA Calculation ---")
print(f"Your total allowed IRA contribution for 2024 (including a ${limit_ira_catchup} catch-up) is: ${total_allowed_ira}")
print(f"You have already contributed: ${user_ira_contribution}")
print(f"Equation: ${total_allowed_ira} (Total Allowed) - ${user_ira_contribution} (Contributed) = ${remaining_ira}")
print(f"Remaining IRA contribution room: ${remaining_ira}\n")

# --- 3. Analyze HSA ---
# The user has already maxed out their HSA for self-only coverage.
remaining_hsa = limit_hsa_self_only - user_hsa_contribution
print("--- HSA & FSA Analysis ---")
print(f"You have contributed ${user_hsa_contribution} to your HSA, which is the maximum for self-only coverage in 2024. Remaining room: ${remaining_hsa}.")
print("FSA contributions do not count toward retirement savings limits.\n")


# --- 4. Calculate Total ---
total_remaining_contribution = remaining_401k + remaining_ira
print("--- Total Remaining Contribution ---")
print("To find your total remaining retirement contribution room, we add the remaining amounts for your 401(k) and IRA.")
print(f"Final Equation: ${remaining_401k} (401k room) + ${remaining_ira} (IRA room) = ${total_remaining_contribution}")
print(f"\nYou can legally contribute a total of ${total_remaining_contribution} more to your retirement accounts in 2024.")
