import datetime

def calculate_deadline():
    """
    Calculates the 5-business-day examination deadline for the bank guarantee.
    """
    presentation_date = datetime.date(2020, 4, 27)
    
    # Vietnamese holidays in that period
    # April 30th (Reunification Day)
    # May 1st (International Labour Day)
    holidays = {
        datetime.date(2020, 4, 30),
        datetime.date(2020, 5, 1)
    }
    
    examination_period_days = 5
    business_days_counted = 0
    
    current_date = presentation_date
    
    print(f"Presentation of documents was on: {presentation_date.strftime('%A, %d %B %Y')}")
    print("--- Calculating 5 business day examination period ---")
    
    # Start counting from the day AFTER presentation
    while business_days_counted < examination_period_days:
        current_date += datetime.timedelta(days=1)
        day_name = current_date.strftime('%A')
        
        # A business day is a weekday and not a holiday.
        # Monday is 0 and Sunday is 6 in weekday()
        if current_date.weekday() < 5 and current_date not in holidays:
            business_days_counted += 1
            print(f"Business Day {business_days_counted}: {current_date.strftime('%d %B %Y')} ({day_name})")
        else:
            reason = "Weekend" if current_date.weekday() >= 5 else "Holiday"
            print(f"Skipping: {current_date.strftime('%d %B %Y')} ({day_name}) - {reason}")
            
    deadline_date = current_date
    print("--- Calculation Complete ---")
    print(f"\nThe examination deadline was the close of business on: {deadline_date.strftime('%A, %d %B %Y')}")
    
    # Final analysis based on the problem description
    print("\nThe issuing bank sent its refusal messages on 06 May 2020, which is within the 5-business-day deadline.")
    print("However, standard banking practice requires a single notice listing all discrepancies.")
    print("The bank sent two separate notices. This procedural failure makes their refusal invalid.")
    print("Therefore, the issuing bank is precluded from claiming non-compliance and must pay.")
    print("The correct action is for the presenting bank to demand that the issuing bank honors the documents.")

calculate_deadline()
