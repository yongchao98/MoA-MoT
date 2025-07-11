# 2024 IRS Contribution Limits
LIMIT_401K_REGULAR = 23000
LIMIT_401K_CATCHUP_AGE_50_PLUS = 7500
LIMIT_IRA_REGULAR = 7000
LIMIT_IRA_CATCHUP_AGE_50_PLUS = 1000
LIMIT_HSA_SELF_ONLY = 4150
LIMIT_FSA_HEALTH = 3200

# User's provided information
user_age_end_of_year = 50
user_contribution_401k = 23000
user_contribution_ira = 5000
user_contribution_hsa = 4150
user_contribution_fsa = 3200

# Since the user turns 50 in 2024, they are eligible for catch-up contributions.
total_limit_401k = LIMIT_401K_REGULAR + LIMIT_401K_CATCHUP_AGE_50_PLUS
total_limit_ira = LIMIT_IRA_REGULAR + LIMIT_IRA_CATCHUP_AGE_50_PLUS

# Calculate remaining contribution room for each retirement account
remaining_401k = total_limit_401k - user_contribution_401k
remaining_ira = total_limit_ira - user_contribution_ira

# The total remaining contribution is the sum of the room in the 401k and IRA.
total_remaining_contribution = remaining_401k + remaining_ira

# Print the breakdown of the calculation
print("--- 2024 Retirement Contribution Analysis ---")
print(f"Your age at the end of 2024 will be {user_age_end_of_year}, making you eligible for catch-up contributions.\n")

print("401k Account:")
print(f"Your total 401k contribution limit is ${total_limit_401k:,} (${LIMIT_401K_REGULAR:,} regular + ${LIMIT_401K_CATCHUP_AGE_50_PLUS:,} catch-up).")
print(f"You have contributed ${user_contribution_401k:,}.")
print(f"Remaining 401k contribution room: ${remaining_401k:,}\n")

print("IRA Account:")
print(f"Your total IRA contribution limit is ${total_limit_ira:,} (${LIMIT_IRA_REGULAR:,} regular + ${LIMIT_IRA_CATCHUP_AGE_50_PLUS:,} catch-up).")
print(f"You have contributed ${user_contribution_ira:,}.")
print(f"Remaining IRA contribution room: ${remaining_ira:,}\n")

print("Other Accounts:")
print(f"Your HSA contribution of ${user_contribution_hsa:,} has reached the annual limit for self-only coverage.")
print(f"Your FSA contribution of ${user_contribution_fsa:,} has reached the annual health FSA limit.\n")

print("--- Total Remaining Contribution ---")
print(f"The total additional amount you can legally contribute to your retirement accounts in 2024 is calculated as follows:")
print(f"${remaining_401k:,} (401k) + ${remaining_ira:,} (IRA) = ${total_remaining_contribution:,}")

<<<10500>>>