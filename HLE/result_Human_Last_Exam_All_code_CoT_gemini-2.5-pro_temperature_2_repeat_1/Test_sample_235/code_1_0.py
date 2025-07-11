def calculate_remaining_retirement_contributions():
    """
    Calculates the remaining legally allowed retirement contributions for 2024
    based on the user's details.
    """

    # User's provided information
    user_age_turns_in_2024 = 50
    user_contrib_401k = 23000
    user_contrib_ira = 5000

    # 2024 IRS Contribution Limits
    limit_401k_standard = 23000
    limit_401k_catchup = 7500
    limit_ira_standard = 7000
    limit_ira_catchup = 1000
    limit_hsa_self = 4150
    limit_fsa = 3200

    print("Retirement Contribution Analysis for 2024\n")

    # --- 401(k) Calculation ---
    print("1. 401(k) Account:")
    # Since the user turns 50 in 2024, they are eligible for the catch-up contribution.
    total_limit_401k = limit_401k_standard + limit_401k_catchup
    remaining_401k = total_limit_401k - user_contrib_401k

    print(f"  Your total 401(k) contribution limit (standard + catch-up) is ${limit_401k_standard:,} + ${limit_401k_catchup:,} = ${total_limit_401k:,}.")
    print(f"  You have already contributed ${user_contrib_401k:,}.")
    print(f"  Remaining 401(k) contribution allowed: ${total_limit_401k:,} - ${user_contrib_401k:,} = ${remaining_401k:,}\n")

    # --- IRA Calculation ---
    print("2. IRA Account (Backdoor Roth):")
    # Since the user turns 50 in 2024, they are eligible for the catch-up contribution.
    total_limit_ira = limit_ira_standard + limit_ira_catchup
    remaining_ira = total_limit_ira - user_contrib_ira

    print(f"  Your total IRA contribution limit (standard + catch-up) is ${limit_ira_standard:,} + ${limit_ira_catchup:,} = ${total_limit_ira:,}.")
    print(f"  You have already contributed ${user_contrib_ira:,}.")
    print(f"  Remaining IRA contribution allowed: ${total_limit_ira:,} - ${user_contrib_ira:,} = ${remaining_ira:,}\n")
    
    # --- Other Accounts ---
    print("3. Other Tax-Advantaged Accounts:")
    print(f"  Your HSA and FSA contributions are already maxed out for 2024 at ${limit_hsa_self:,} and ${limit_fsa:,} respectively.\n")

    # --- Final Calculation ---
    total_remaining = remaining_401k + remaining_ira
    print("--- Summary ---")
    print(f"Total additional money you are legally allowed to contribute to your retirement accounts in 2024 is:")
    print(f"Remaining 401(k) + Remaining IRA = Total")
    print(f"${remaining_401k:,} + ${remaining_ira:,} = ${total_remaining:,}")

    return total_remaining

final_answer = calculate_remaining_retirement_contributions()
print(f"<<<${final_answer:,.0f}>>>")