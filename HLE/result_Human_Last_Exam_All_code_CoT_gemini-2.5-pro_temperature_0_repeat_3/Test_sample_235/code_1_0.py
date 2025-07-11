# 2024 IRS Contribution Limits
MAX_401K_REGULAR = 23000
MAX_401K_CATCHUP = 7500
MAX_IRA_REGULAR = 7000
MAX_IRA_CATCHUP = 1000

# User's current contributions
user_401k_contrib = 23000
user_ira_contrib = 5000

# Since you turn 50 in 2024, you are eligible for catch-up contributions.
total_401k_limit = MAX_401K_REGULAR + MAX_401K_CATCHUP
total_ira_limit = MAX_IRA_REGULAR + MAX_IRA_CATCHUP

# Calculate remaining contribution room
remaining_401k = total_401k_limit - user_401k_contrib
remaining_ira = total_ira_limit - user_ira_contrib

# Calculate total remaining contributions
total_remaining = remaining_401k + remaining_ira

# Print the step-by-step calculation
print("--- 401(k) Calculation ---")
print(f"Your total 401(k) contribution limit (including age 50+ catch-up) is: ${MAX_401K_REGULAR} + ${MAX_401K_CATCHUP} = ${total_401k_limit}")
print(f"You have already contributed: ${user_401k_contrib}")
print(f"Remaining 401(k) contribution allowed: ${total_401k_limit} - ${user_401k_contrib} = ${remaining_401k}\n")

print("--- IRA Calculation ---")
print(f"Your total IRA contribution limit (including age 50+ catch-up) is: ${MAX_IRA_REGULAR} + ${MAX_IRA_CATCHUP} = ${total_ira_limit}")
print(f"You have already contributed: ${user_ira_contrib}")
print(f"Remaining IRA contribution allowed: ${total_ira_limit} - ${user_ira_contrib} = ${remaining_ira}\n")

print("--- Total Additional Contribution ---")
print("The total additional amount you can contribute to your retirement accounts in 2024 is:")
print(f"${remaining_401k} (401k) + ${remaining_ira} (IRA) = ${total_remaining}")