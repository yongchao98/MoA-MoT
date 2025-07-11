def calculate_remaining_retirement_contributions():
    """
    Calculates and explains the remaining retirement contributions for a 49-year-old
    turning 50 in the 2024 tax year.
    """
    # --- 2024 IRS Limits ---
    limit_401k_standard = 23000
    limit_401k_catchup_age_50_plus = 7500
    limit_ira_standard = 7000
    limit_ira_catchup_age_50_plus = 1000
    limit_hsa_self_only = 4150
    limit_fsa_health = 3200

    # --- User's Contributions ---
    user_contrib_401k = 23000
    user_contrib_ira = 5000
    user_contrib_hsa = 4150
    user_contrib_fsa = 3200

    # --- Explanation ---
    print("Based on the 2024 tax year rules, here is the breakdown of your remaining contribution room:")
    print("-" * 70)

    # Step 1: 401(k) Calculation
    # The user turns 50 in 2024, making them eligible for catch-up contributions.
    total_401k_limit = limit_401k_standard + limit_401k_catchup_age_50_plus
    remaining_401k = total_401k_limit - user_contrib_401k
    print(f"1. 401(k) Account:")
    print(f"   - The standard employee contribution limit for 2024 is ${limit_401k_standard:,}.")
    print(f"   - Because you turn 50 in 2024, you are eligible for an additional catch-up contribution of ${limit_401k_catchup_age_50_plus:,}.")
    print(f"   - Your total 401(k) contribution limit is ${total_401k_limit:,}.")
    print(f"   - You have already contributed ${user_contrib_401k:,}.")
    print(f"   - Amount you can still contribute to your 401(k): ${remaining_401k:,}")
    print("\n   Note: Employer contributions do not count against your personal contribution limit.\n")

    # Step 2: IRA Calculation
    total_ira_limit = limit_ira_standard + limit_ira_catchup_age_50_plus
    remaining_ira = total_ira_limit - user_contrib_ira
    print(f"2. IRA Account (Backdoor Roth):")
    print(f"   - The standard IRA contribution limit for 2024 is ${limit_ira_standard:,}.")
    print(f"   - Your age also makes you eligible for the IRA catch-up contribution of ${limit_ira_catchup_age_50_plus:,}.")
    print(f"   - Your total IRA contribution limit is ${total_ira_limit:,}.")
    print(f"   - You have already contributed ${user_contrib_ira:,}.")
    print(f"   - Amount you can still contribute to your IRA: ${remaining_ira:,}\n")
    
    # Step 3: Non-Retirement Accounts
    remaining_hsa = limit_hsa_self_only - user_contrib_hsa
    remaining_fsa = limit_fsa_health - user_contrib_fsa
    print("3. Other Accounts (Not included in retirement total):")
    print(f"   - Health Savings Account (HSA): You have contributed ${user_contrib_hsa:,}, which is the maximum for self-only coverage in 2024. Remaining: ${remaining_hsa:,}.")
    print(f"   - Flexible Spending Account (FSA): You have contributed ${user_contrib_fsa:,}, which is the maximum for 2024. Remaining: ${remaining_fsa:,}.\n")

    # Final Calculation
    total_remaining = remaining_401k + remaining_ira
    print("-" * 70)
    print("Final Calculation for Total Remaining Retirement Contributions:\n")
    print("Remaining 401k Contribution + Remaining IRA Contribution = Total")
    # Using format specifiers to ensure proper output for the equation
    print(f"${remaining_401k:,} + ${remaining_ira:,} = ${total_remaining:,}")


calculate_remaining_retirement_contributions()