def calculate_remaining_contributions():
    """
    Calculates and explains the remaining retirement contributions for 2024
    based on the user's provided information.
    """
    # 2024 Contribution Limits
    LIMIT_401K_STANDARD = 23000
    LIMIT_401K_CATCHUP = 7500
    LIMIT_IRA_STANDARD = 7000
    LIMIT_IRA_CATCHUP = 1000
    LIMIT_HSA_SELF = 4150

    # User's Contributions
    contrib_401k = 23000
    contrib_ira = 5000
    contrib_hsa = 4150

    # The user turns 50 in 2024, making them eligible for catch-up contributions.
    total_401k_limit = LIMIT_401K_STANDARD + LIMIT_401K_CATCHUP
    total_ira_limit = LIMIT_IRA_STANDARD + LIMIT_IRA_CATCHUP
    total_hsa_limit = LIMIT_HSA_SELF # HSA catch-up starts at age 55

    # Calculate remaining contribution room for each account type
    remaining_401k = total_401k_limit - contrib_401k
    remaining_ira = total_ira_limit - contrib_ira
    remaining_hsa = total_hsa_limit - contrib_hsa
    
    total_remaining = remaining_401k + remaining_ira + remaining_hsa

    print("Here is a breakdown of your remaining tax-advantaged account contributions for 2024.\n")
    print("Note: Since you turn 50 in 2024, you are eligible for catch-up contributions for both your 401k and IRA.")

    print("\n--- 401k ---")
    print(f"Your total 401k contribution limit is ${total_401k_limit:,} (${LIMIT_401K_STANDARD:,} standard + ${LIMIT_401K_CATCHUP:,} catch-up).")
    print(f"You have contributed ${contrib_401k:,}.")
    print(f"Calculation: ${total_401k_limit:,} (Limit) - ${contrib_401k:,} (Contributed) = ${remaining_401k:,} remaining.")

    print("\n--- IRA (Backdoor Roth) ---")
    print(f"Your total IRA contribution limit is ${total_ira_limit:,} (${LIMIT_IRA_STANDARD:,} standard + ${LIMIT_IRA_CATCHUP:,} catch-up).")
    print(f"You have contributed ${contrib_ira:,}.")
    print(f"Calculation: ${total_ira_limit:,} (Limit) - ${contrib_ira:,} (Contributed) = ${remaining_ira:,} remaining.")
    
    print("\n--- HSA (Health Savings Account) ---")
    print(f"Your HSA contribution limit for self-only coverage is ${total_hsa_limit:,}.")
    print(f"You have contributed ${contrib_hsa:,}.")
    print(f"Calculation: ${total_hsa_limit:,} (Limit) - ${contrib_hsa:,} (Contributed) = ${remaining_hsa:,} remaining.")
    print("(Note: HSA catch-up contributions begin at age 55.)")

    print("\n--- FSA & Employer Match ---")
    print("Your FSA and employer's 401k match do not count against your personal contribution limits and are excluded from this calculation.")

    print("\n-------------------------------------------")
    print("TOTAL REMAINING CONTRIBUTIONS FOR 2024")
    print("-------------------------------------------")
    print(f"Remaining 401k (${remaining_401k:,}) + Remaining IRA (${remaining_ira:,}) + Remaining HSA (${remaining_hsa:,}) = ${total_remaining:,}")
    print("-------------------------------------------")

calculate_remaining_contributions()