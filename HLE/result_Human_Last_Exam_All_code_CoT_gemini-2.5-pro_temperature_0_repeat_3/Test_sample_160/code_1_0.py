import datetime

def solve_trade_finance_case():
    """
    Analyzes the trade finance case based on URDG 758 rules.
    Calculates the examination period and determines the correct action.
    """
    presentation_date = datetime.date(2020, 4, 27)
    holidays = [
        datetime.date(2020, 4, 30),  # Vietnamese holiday
        datetime.date(2020, 5, 1)   # International Labour Day
    ]
    
    print("Analyzing the examination period based on URDG 758 Article 24...")
    print(f"Presentation of documents was on: {presentation_date.strftime('%A, %d %B %Y')}")
    print("The issuing bank has a maximum of five business days following the day of presentation to send a notice of refusal.")
    print("-" * 60)

    business_days_passed = 0
    current_date = presentation_date

    while business_days_passed < 5:
        current_date += datetime.timedelta(days=1)
        day_of_week = current_date.weekday() # Monday is 0, Sunday is 6

        # Check if it's a weekend (Saturday or Sunday)
        if day_of_week >= 5:
            status = "Weekend"
        # Check if it's a holiday
        elif current_date in holidays:
            status = "Holiday"
        # Otherwise, it's a business day
        else:
            business_days_passed += 1
            status = f"Business Day #{business_days_passed}"
        
        print(f"Date: {current_date.strftime('%A, %d %B %Y')} - Status: {status}")

    deadline_date = current_date
    print("-" * 60)
    print(f"The deadline for the issuing bank to send a refusal notice was the close of business on: {deadline_date.strftime('%A, %d %B %Y')}")
    
    print("\nAnalyzing the Notice of Refusal:")
    print("1. The issuing bank sent a first refusal notice on 06 May 2020 at 09:00 am.")
    print("2. The issuing bank sent a second refusal notice on 06 May 2020 at 03:00 pm with a new discrepancy.")
    
    print("\nApplying URDG 758 Rules:")
    print("URDG 758 Article 24(c) requires a 'single notice of refusal' that lists 'all discrepancies'.")
    print("By sending two separate notices, the issuing bank failed to comply with this rule.")
    print("URDG 758 Article 24(f) states that if a bank fails to comply, it is 'precluded from claiming that the demand...[is] not in conformity'.")

    print("\nConclusion:")
    print("The issuing bank has lost its right to refuse the demand due to its procedural error.")
    print("Therefore, the presenting bank should insist on payment.")

solve_trade_finance_case()
<<<A>>>