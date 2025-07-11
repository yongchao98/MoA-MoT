import datetime

def solve_bank_guarantee_case():
    """
    Analyzes a bank guarantee case based on URDG 760 rules,
    calculates the examination period, and provides a recommendation.
    """
    # --- Define Key Dates and Parameters ---
    presentation_date = datetime.date(2020, 4, 27)
    first_refusal_date = datetime.date(2020, 5, 6)
    
    # Vietnamese holidays in that period
    # 30 April: Reunification Day
    # 01 May: Labour Day
    holidays = [
        datetime.date(2020, 4, 30),
        datetime.date(2020, 5, 1)
    ]
    
    urdeg_examination_days = 5

    # --- Step 1: Calculate the Examination Period Deadline ---
    print("--- Step 1: Calculating the Examination Period Deadline (URDG 760 Article 20) ---")
    print(f"A demand was presented on: {presentation_date.strftime('%d %b %Y, %A')}")
    print(f"The issuing bank has {urdeg_examination_days} business days from the following day to examine the demand.")
    
    business_days_counted = 0
    current_date = presentation_date
    examination_days_list = []

    # Loop to find the 5 business days
    while business_days_counted < urdeg_examination_days:
        current_date += datetime.timedelta(days=1)
        # weekday() -> Monday is 0 and Sunday is 6
        is_weekend = current_date.weekday() >= 5
        is_holiday = current_date in holidays
        
        day_info = f"{current_date.strftime('%d %b %Y, %A')}"

        if not is_weekend and not is_holiday:
            business_days_counted += 1
            examination_days_list.append(f"Business Day {business_days_counted}: {day_info}")
        elif is_weekend:
            # This is just for verbose output, not added to the list
            pass
        elif is_holiday:
            # This is just for verbose output, not added to the list
            pass
            
    deadline_date = current_date

    print("\nThe 5 business days for examination are:")
    for day in examination_days_list:
        print(f"- {day}")

    print(f"\nThe deadline for the issuing bank to send a refusal was the close of business on: {deadline_date.strftime('%d %b %Y')}")
    print(f"The issuing bank sent its first refusal on: {first_refusal_date.strftime('%d %b %Y')}")

    if first_refusal_date <= deadline_date:
        print("Conclusion: The refusal was sent *within* the allowed time.")
    else:
        print("Conclusion: The refusal was sent *after* the allowed time.")

    # --- Step 2: Analyze the Content of the Refusal (URDG 760 Article 24) ---
    print("\n--- Step 2: Analyzing the Content of the Refusal (URDG 760 Article 24) ---")
    print("URDG Article 24(c) requires that a *single notice* of refusal must state *all discrepancies*.")
    print("\nThe issuing bank sent two separate notices on 06 May 2020:")
    print("1. Notice at 09:00 am: Discrepancy of wrong demand amount.")
    print("2. Notice at 03:00 pm: Discrepancy of missing beneficiary's signature.")

    print("\nURDG Article 24(d) states that if a bank fails to act in accordance with paragraph (c), it is precluded (prevented) from claiming the demand is non-complying.")
    print("By sending a first notice that did not list all discrepancies, the bank failed to follow the rules.")
    print("Therefore, the issuing bank can only rely on the discrepancy mentioned in the first notice.")
    
    print("\n--- Final Conclusion ---")
    print("Valid Discrepancy: 'Wrong demand amount'.")
    print("Invalid Discrepancy: 'Missing beneficiary's signature' (because it was raised in a second, separate notice).")
    print("\nThe presenting bank should ignore the second discrepancy and advise the beneficiary about the first, valid discrepancy so it can be corrected.")

solve_bank_guarantee_case()