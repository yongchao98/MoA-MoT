def calculate_remaining_retirement_contributions():
    """
    Calculates the remaining retirement contribution allowance for a 49-year-old
    turning 50 in the 2024 tax year.
    """
    # IRS limits for 2024
    REGULAR_401K_LIMIT = 23000
    CATCH_UP_401K = 7500
    REGULAR_IRA_LIMIT = 7000
    CATCH_UP_IRA = 1000

    # User's provided contributions
    user_401k_contribution = 23000
    user_ira_contribution = 5000

    # As the user will turn 50 in 2024, they are eligible for catch-up contributions.
    total_401k_limit = REGULAR_401K_LIMIT + CATCH_UP_401K
    total_ira_limit = REGULAR_IRA_LIMIT + CATCH_UP_IRA

    # Calculate remaining contribution room for each account
    remaining_401k = total_401k_limit - user_401k_contribution
    remaining_ira = total_ira_limit - user_ira_contribution

    # Calculate total remaining contribution
    total_remaining_contribution = remaining_401k + remaining_ira

    print("Based on the information provided for the 2024 tax year, here is your remaining retirement contribution allowance.")
    print("-" * 70)
    print("Note: Since you turn age 50 in 2024, you are eligible for 'catch-up' contributions.")
    print("Your HSA ($4,150) and FSA ($3,200) are already maxed out for the year.")
    print("-" * 70)

    print("\n401(k) Contribution Analysis:")
    print(f"Your total 401(k) contribution limit (standard + catch-up) is: ${REGULAR_401K_LIMIT} + ${CATCH_UP_401K} = ${total_401k_limit}")
    print(f"You have contributed: ${user_401k_contribution}")
    print(f"Remaining 401(k) contribution room: ${total_401k_limit} - ${user_401k_contribution} = ${remaining_401k}")
    
    print("\nIRA Contribution Analysis:")
    print(f"Your total IRA contribution limit (standard + catch-up) is: ${REGULAR_IRA_LIMIT} + ${CATCH_UP_IRA} = ${total_ira_limit}")
    print(f"You have contributed: ${user_ira_contribution}")
    print(f"Remaining IRA contribution room: ${total_ira_limit} - ${user_ira_contribution} = ${remaining_ira}")
    print("-" * 70)
    
    print("\nTotal Remaining Retirement Contributions for 2024:")
    print(f"The total additional amount you can contribute is calculated as:")
    print(f"${remaining_401k} (to your 401k) + ${remaining_ira} (to your IRA) = ${total_remaining_contribution}")

calculate_remaining_retirement_contributions()
print(f'<<<{10500}>>>')