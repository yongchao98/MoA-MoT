import datetime
from datetime import timedelta

def solve_guarantee_case():
    """
    Analyzes the bank guarantee case by calculating the examination period deadline
    and interpreting the rules of URDG 758.
    """
    presentation_date = datetime.date(2020, 4, 27)
    
    # Holidays in Vietnam for that period
    # April 30 (Reunification Day), May 1 (Labour Day)
    holidays = {
        datetime.date(2020, 4, 30),
        datetime.date(2020, 5, 1)
    }

    print("Step 1: Calculating the Examination Period Deadline based on URDG 758 Article 20.")
    print(f"Presentation was made on: {presentation_date.strftime('%A, %d %B %Y')}")
    print("The examination period is a maximum of five business days following the day of presentation.\n")

    business_days_counted = 0
    current_date = presentation_date
    
    # Loop to find the 5 business days
    while business_days_counted < 5:
        current_date += timedelta(days=1)
        # weekday() returns 0 for Monday and 6 for Sunday. We need Monday-Friday.
        is_weekday = current_date.weekday() < 5
        is_holiday = current_date in holidays
        
        if is_weekday and not is_holiday:
            business_days_counted += 1
            print(f"Business Day {business_days_counted}: {current_date.strftime('%A, %d %B %Y')}")
        else:
            reason = "Weekend" if not is_weekday else "Holiday"
            print(f"Skipping           : {current_date.strftime('%A, %d %B %Y')} ({reason})")

    deadline_date = current_date
    print(f"\nStep 2: Determine the deadline.")
    print(f"The deadline for the issuing bank to send a refusal notice was the close of business on: {deadline_date.strftime('%A, %d %B %Y')}.")

    refusal_date_1 = datetime.datetime(2020, 5, 6, 9, 0)
    refusal_date_2 = datetime.datetime(2020, 5, 6, 15, 0)
    
    print(f"\nThe bank sent its first refusal message on {refusal_date_1.strftime('%d %B %Y, at %H:%M')}.")
    print("This was within the 5-business-day time limit.")
    print("\nStep 3: Analyze the content and method of refusal based on URDG 758 Article 24.")
    print("URDG 758 Article 24(c) requires a notice of rejection to state ALL discrepancies in a single notice.")
    print(f"The issuing bank sent a first message at 09:00 am citing one discrepancy (wrong amount).")
    print(f"It then sent a SECOND message at 03:00 pm on the same day citing another discrepancy (missing signature).")
    print("\nStep 4: Conclusion based on the rules.")
    print("By sending two separate notices, the issuing bank violated URDG 758 Article 24(c).")
    print("According to URDG 758 Article 24(d), this failure means the bank is 'precluded from claiming that the demand... is not complying.'")
    print("Therefore, the issuing bank has lost its right to refuse payment and must honor the demand.")
    print("\nThe correct action for the presenting bank is to point out this procedural failure and demand payment.")

solve_guarantee_case()
<<<A>>>