import datetime

def solve_lc_deadline():
    """
    Calculates the deadline for a bank to send a refusal message for a Letter of Credit presentation.
    """
    # Step 1: Define initial parameters
    presentation_datetime_str = "2020-03-30 17:01"
    presentation_datetime = datetime.datetime.strptime(presentation_datetime_str, "%Y-%m-%d %H:%M")
    bank_closing_time = datetime.time(17, 0)
    # 5 = Saturday, 6 = Sunday
    non_banking_days = [5, 6]  
    examination_period_days = 5

    print(f"Initial Presentation: {presentation_datetime.strftime('%A, %d %B %Y, at %I:%M %p')}")
    print(f"Bank's Closing Time: {bank_closing_time.strftime('%I:%M %p')}")
    print("Rule: If presentation is after closing time, it's considered received on the next banking day.\n")

    # Step 2: Determine the effective date of presentation
    effective_presentation_date = presentation_datetime.date()
    if presentation_datetime.time() > bank_closing_time:
        print(f"Presentation time ({presentation_datetime.strftime('%I:%M %p')}) is after the closing time.")
        
        # Find the next banking day
        next_day = effective_presentation_date + datetime.timedelta(days=1)
        while next_day.weekday() in non_banking_days:
            next_day += datetime.timedelta(days=1)
        effective_presentation_date = next_day
        print(f"Therefore, the effective day of presentation is the next banking day: {effective_presentation_date.strftime('%A, %d %B %Y')}")
    else:
        print(f"Presentation time is within business hours. Effective date is {effective_presentation_date.strftime('%A, %d %B %Y')}")

    print("\n-------------------------------------------------------------")
    print(f"Rule: The bank has a maximum of {examination_period_days} banking days FOLLOWING the day of presentation to send a refusal.")
    print(f"Calculation starts from the day after {effective_presentation_date.strftime('%A, %d %B %Y')}.\n")

    # Step 3: Count the 5 banking days
    current_date = effective_presentation_date
    days_counted = 0
    
    while days_counted < examination_period_days:
        current_date += datetime.timedelta(days=1)
        if current_date.weekday() not in non_banking_days:
            days_counted += 1
            print(f"Day {days_counted}: {current_date.strftime('%A, %d %B %Y')} (Banking Day)")
        else:
            print(f"      {current_date.strftime('%A, %d %B %Y')} (Non-Banking Day - Skipped)")
    
    deadline_date = current_date
    print("\n-------------------------------------------------------------")
    print(f"The fifth and final banking day is {deadline_date.strftime('%A, %d %B %Y')}.")
    print(f"Therefore, the bank should send the refusal message at the latest on {deadline_date.strftime('%d %B %Y')}.")


solve_lc_deadline()
<<<D>>>