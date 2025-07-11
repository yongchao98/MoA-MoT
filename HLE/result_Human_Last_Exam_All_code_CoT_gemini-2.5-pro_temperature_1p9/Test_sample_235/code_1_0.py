def calculate_remaining_contributions():
    """
    Calculates the remaining legal retirement contributions for a 49-year-old in the 2024 tax year.
    """
    # User's current contributions
    user_401k_contribution = 23000
    user_ira_contribution = 5000
    user_hsa_contribution = 4150
    user_fsa_contribution = 3200

    # 2024 IRS Limits
    limit_401k_standard = 23000
    limit_401k_catchup = 7500  # For age 50+
    limit_ira_standard = 7000
    limit_ira_catchup = 1000  # For age 50+
    limit_hsa_self_only = 4150
    limit_fsa_health = 3200

    # User is turning 50 in 2024, so they are eligible for catch-up contributions for the full year.
    total_allowed_401k = limit_401k_standard + limit_401k_catchup
    total_allowed_ira = limit_ira_standard + limit_ira_catchup

    # Calculate remaining contribution room
    remaining_401k = total_allowed_401k - user_401k_contribution
    remaining_ira = total_allowed_ira - user_ira_contribution
    total_remaining_retirement = remaining_401k + remaining_ira

    print("--- 2024 Retirement Contribution Analysis ---")
    print("\nSince you turn 50 in 2024, you are eligible for catch-up contributions for the entire year.")

    # 401(k) Calculation
    print("\n--- 401(k) Account ---")
    print(f"Standard employee limit: ${limit_401k_standard}")
    print(f"Age 50+ catch-up limit: ${limit_401k_catchup}")
    print(f"Your total allowable 401(k) contribution is ${limit_401k_standard} + ${limit_401k_catchup} = ${total_allowed_401k}")
    print(f"You have already contributed: ${user_401k_contribution}")
    print(f"Remaining 401(k) contribution for 2024: ${remaining_401k}")

    # IRA Calculation
    print("\n--- IRA (Backdoor Roth) Account ---")
    print(f"Standard IRA limit: ${limit_ira_standard}")
    print(f"Age 50+ catch-up limit: ${limit_ira_catchup}")
    print(f"Your total allowable IRA contribution is ${limit_ira_standard} + ${limit_ira_catchup} = ${total_allowed_ira}")
    print(f"You have already contributed: ${user_ira_contribution}")
    print(f"Remaining IRA contribution for 2024: ${remaining_ira}")

    # Other Accounts
    print("\n--- Other Accounts ---")
    print(f"Your HSA contribution of ${user_hsa_contribution} has reached the annual maximum of ${limit_hsa_self_only}.")
    print(f"Your FSA contribution of ${user_fsa_contribution} has reached the annual maximum of ${limit_fsa_health}.")
    
    # Final Summary
    print("\n--- Total Remaining Contribution ---")
    print("The total additional money you can legally contribute to your retirement accounts in 2024 is:")
    print(f"Remaining 401(k) Room + Remaining IRA Room")
    print(f"${remaining_401k} + ${remaining_ira} = ${total_remaining_retirement}")

calculate_remaining_contributions()