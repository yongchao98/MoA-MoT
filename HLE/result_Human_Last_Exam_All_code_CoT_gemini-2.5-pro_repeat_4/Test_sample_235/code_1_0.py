def calculate_remaining_retirement_contributions():
    """
    Calculates the remaining amount a user can legally contribute to retirement accounts in 2024.
    """
    # User's information
    user_age = 49  # Turning 50 in 2024, making them eligible for catch-up contributions
    user_401k_contribution = 23000
    user_ira_contribution = 5000
    user_hsa_contribution = 4150
    
    # 2024 IRS limits
    limit_401k_standard = 23000
    limit_401k_catch_up = 7500
    limit_ira_standard = 7000
    limit_ira_catch_up = 1000
    limit_hsa_single = 4150

    # Since the user turns 50 in 2024, they are eligible for catch-up contributions for the entire year.
    total_401k_limit = limit_401k_standard + limit_401k_catch_up
    total_ira_limit = limit_ira_standard + limit_ira_catch_up

    # Calculate remaining contribution room for each retirement account
    remaining_401k = total_401k_limit - user_401k_contribution
    remaining_ira = total_ira_limit - user_ira_contribution
    
    # The HSA is already maxed out for 2024 ($4,150 limit - $4,150 contributed = $0).
    # The FSA is not a retirement account.
    
    # Calculate total remaining for retirement accounts
    total_remaining = remaining_401k + remaining_ira
    
    print("Calculating your remaining retirement contributions for 2024:\n")
    
    print("401k Calculation:")
    print(f"Your total 401k limit (including age 50+ catch-up): ${limit_401k_standard} + ${limit_401k_catch_up} = ${total_401k_limit}")
    print(f"Remaining 401k Contribution: ${total_401k_limit} - ${user_401k_contribution} = ${remaining_401k}\n")
    
    print("IRA Calculation:")
    print(f"Your total IRA limit (including age 50+ catch-up): ${limit_ira_standard} + ${limit_ira_catch_up} = ${total_ira_limit}")
    print(f"Remaining IRA Contribution: ${total_ira_limit} - ${user_ira_contribution} = ${remaining_ira}\n")

    print("Total Remaining Retirement Contributions Calculation:")
    print(f"Remaining 401k (${remaining_401k}) + Remaining IRA (${remaining_ira})")
    print("========================================")
    print(f"Total: ${total_remaining}")

calculate_remaining_retirement_contributions()
<<<10500>>>