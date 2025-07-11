import datetime

def calculate_refusal_deadline():
    """
    Calculates the deadline for an issuing bank to refuse a demand under URDG 758.
    
    The script identifies the 5 business days following the presentation,
    excluding weekends and public holidays, to determine the final deadline.
    """
    
    # --- Step 1: Define key dates and parameters ---
    presentation_date = datetime.date(2020, 4, 27)
    
    # Vietnamese holidays in the period.
    # April 30th is Reunification Day. May 1st is International Labor Day.
    holidays = [
        datetime.date(2020, 4, 30),
        datetime.date(2020, 5, 1)
    ]
    
    examination_period_days = 5
    business_days_counted = 0
    current_date = presentation_date

    print(f"Analysis based on URDG 758 Article 24")
    print("-----------------------------------------")
    print(f"Presentation Date: {presentation_date.strftime('%A, %d %B %Y')}")
    print(f"Examination period is 5 business days following the day of presentation.")
    print("-----------------------------------------\n")
    print("Calculating the deadline:\n")
    
    # --- Step 2: Loop to find the 5 business days ---
    while business_days_counted < examination_period_days:
        current_date += datetime.timedelta(days=1)
        day_name = current_date.strftime('%A')
        
        # Saturday is 5, Sunday is 6 in weekday() method
        if current_date.weekday() in [5, 6]:
            print(f"{current_date.strftime('%d %B %Y')} ({day_name}) - Skipped (Weekend)")
            continue
        
        if current_date in holidays:
            print(f"{current_date.strftime('%d %B %Y')} ({day_name}) - Skipped (Public Holiday)")
            continue
        
        business_days_counted += 1
        print(f"Business Day {business_days_counted}: {current_date.strftime('%d %B %Y')} ({day_name})")

    # --- Step 3: State the conclusion ---
    deadline_date = current_date
    refusal_date = datetime.date(2020, 5, 6)
    
    print("\n-----------------------------------------")
    print(f"Conclusion of Calculation")
    print("-----------------------------------------")
    print(f"The deadline for the issuing bank to send a refusal was the end of the business day on: {deadline_date.strftime('%A, %d %B %Y')}")
    print(f"The bank sent its first refusal message on: {refusal_date.strftime('%A, %d %B %Y')}")
    print("\nFinal Assessment:")
    if refusal_date > deadline_date:
        print("The bank sent its refusal AFTER the 5-business-day deadline.")
        print("According to URDG 758 Article 24(f), the issuing bank is precluded from claiming the presentation is non-complying.")
        print("Therefore, the issuing bank has lost its right to refuse and must honor the demand.")
    else:
        print("The bank sent its refusal within the 5-business-day deadline.")

if __name__ == '__main__':
    calculate_refusal_deadline()