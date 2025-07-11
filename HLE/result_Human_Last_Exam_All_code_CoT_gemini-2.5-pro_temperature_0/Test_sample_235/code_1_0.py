def calculate_remaining_retirement_contributions():
    """
    Calculates the remaining legal retirement contributions for a user in 2024.
    """
    # User's information
    user_age = 49  # Turns 50 in 2024, making them eligible for catch-up contributions
    user_401k_contribution = 23000
    user_ira_contribution = 5000
    user_hsa_contribution = 4150

    # 2024 IRS Limits
    MAX_401K_REGULAR = 23000
    MAX_401K_CATCHUP = 7500
    MAX_IRA_REGULAR = 7000
    MAX_IRA_CATCHUP = 1000
    MAX_HSA_SINGLE = 4150
    CATCHUP_ELIGIBILITY_AGE = 50

    # --- 401k Calculation ---
    # User is eligible for catch-up since they turn 50 in the tax year.
    total_401k_limit = MAX_401K_REGULAR + MAX_401K_CATCHUP
    remaining_401k = total_401k_limit - user_401k_contribution

    # --- IRA Calculation ---
    # User is also eligible for the IRA catch-up.
    total_ira_limit = MAX_IRA_REGULAR + MAX_IRA_CATCHUP
    remaining_ira = total_ira_limit - user_ira_contribution

    # --- HSA Calculation ---
    # HSA catch-up starts at age 55.
    total_hsa_limit = MAX_HSA_SINGLE
    remaining_hsa = total_hsa_limit - user_hsa_contribution
    # Ensure it doesn't go below zero if over-contributed
    remaining_hsa = max(0, remaining_hsa)

    # --- Total Calculation ---
    total_remaining = remaining_401k + remaining_ira + remaining_hsa

    # --- Print the explanation ---
    print("Based on your age and 2024 tax laws, here is your remaining retirement contribution room:\n")

    print("401k Account:")
    print(f"  - Your total 401k contribution limit (including a ${MAX_401K_CATCHUP:,} catch-up): ${total_401k_limit:,}")
    print(f"  - You have already contributed: ${user_401k_contribution:,}")
    print(f"  - Remaining 401k contribution room: ${remaining_401k:,}\n")

    print("IRA Account (Backdoor Roth):")
    print(f"  - Your total IRA contribution limit (including a ${MAX_IRA_CATCHUP:,} catch-up): ${total_ira_limit:,}")
    print(f"  - You have already contributed: ${user_ira_contribution:,}")
    print(f"  - Remaining IRA contribution room: ${remaining_ira:,}\n")

    print("HSA Account:")
    print(f"  - Your total HSA contribution limit is: ${total_hsa_limit:,}")
    print(f"  - You have already contributed: ${user_hsa_contribution:,}")
    print(f"  - Remaining HSA contribution room: ${remaining_hsa:,}\n")

    print("Note: Your FSA contributions do not count towards retirement savings limits.")
    print("-" * 40)
    print(f"Total additional amount you can contribute to these accounts in 2024: ${total_remaining:,}")
    print("Final Equation:")
    print(f"${remaining_401k:,} (401k) + ${remaining_ira:,} (IRA) + ${remaining_hsa:,} (HSA) = ${total_remaining:,}")

calculate_remaining_retirement_contributions()
<<<10500>>>