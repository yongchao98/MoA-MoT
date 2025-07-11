import datetime

def analyze_guarantee_case():
    """
    Analyzes the bank guarantee scenario based on URDG 758 rules.
    """
    # 1. Define dates and constants from the problem description
    presentation_date_str = "2020-04-27"
    rejection_date_str = "2020-05-06"
    holidays_str = [
        "2020-04-30", # Vietnamese holiday (Reunification Day)
        "2020-05-01"  # Vietnamese holiday (International Labour Day)
    ]
    
    # Convert strings to date objects for calculation
    presentation_date = datetime.datetime.strptime(presentation_date_str, "%Y-%m-%d").date()
    rejection_date = datetime.datetime.strptime(rejection_date_str, "%Y-%m-%d").date()
    holidays = [datetime.datetime.strptime(d, "%Y-%m-%d").date() for d in holidays_str]
    
    examination_period_days = 5

    # --- Step 1: Calculate the Examination Deadline ---
    print("--- Analysis of the Issuing Bank's Examination Period ---")
    print(f"Presentation of documents was made on: {presentation_date.strftime('%A, %d %B %Y')}")
    print(f"Under URDG 758, the issuing bank has {examination_period_days} business days following the day of presentation to examine the documents.")
    print("Calculating the deadline by checking each day:")

    # The examination period starts on the business day *following* the day of presentation
    current_date = presentation_date + datetime.timedelta(days=1)
    business_days_counted = 0
    examination_deadline = None
    
    while business_days_counted < examination_period_days:
        # weekday() returns Monday as 0 and Sunday as 6
        is_weekend = current_date.weekday() >= 5 # True for Saturday (5) or Sunday (6)
        is_holiday = current_date in holidays
        
        day_status = ""
        if not is_weekend and not is_holiday:
            business_days_counted += 1
            day_status = f"(Business Day {business_days_counted})"
            if business_days_counted == examination_period_days:
                examination_deadline = current_date
        elif is_weekend:
            day_status = "(Weekend - Non-Banking Day)"
        else: # is_holiday
            day_status = "(Holiday - Non-Banking Day)"
            
        print(f" - {current_date.strftime('%d %B %Y, %A')}: {day_status}")
        current_date += datetime.timedelta(days=1)
        
    print(f"\nResult: The deadline for refusal was the close of business on {examination_deadline.strftime('%A, %d %B %Y')}.")
    print(f"The bank sent its first refusal notice on {rejection_date.strftime('%A, %d %B %Y')}, which is within this time limit.")

    # --- Step 2: Analyze the Content and Method of Refusal ---
    print("\n--- Analysis of the Refusal Notice(s) ---")
    print("URDG 758 Article 24(c)(i) is critical here. It requires a refusal notice to state 'each discrepancy' for which the bank rejects the demand.")
    print("This means all discrepancies must be listed in a single, comprehensive notice.")
    print("\nThe issuing bank's actions:")
    print(" 1. At 09:00 am on 06 May, it sent a notice for one discrepancy: 'wrong demand amount'.")
    print(" 2. At 03:00 pm on 06 May, it sent a second notice for another discrepancy: 'missing beneficiary’s signature'.")
    print("\nThis is a 'piecemeal rejection'. By sending the first notice with only one discrepancy, the bank effectively waived all other discrepancies.")
    print("According to URDG 758 Article 24(d), a bank that fails to act in accordance with these rules is 'precluded' from claiming the presentation is non-compliant.")

    # --- Step 3: Conclusion and Recommended Action ---
    print("\n--- Final Conclusion ---")
    print("Because the issuing bank issued a defective, piecemeal rejection, it has lost its right to refuse the demand.")
    print("Therefore, the issuing bank is obligated to honor the presentation, despite the original discrepancies.")
    print("The correct course of action for the presenting bank is to point out the defective refusal and demand payment.")

# Run the analysis
analyze_guarantee_case()