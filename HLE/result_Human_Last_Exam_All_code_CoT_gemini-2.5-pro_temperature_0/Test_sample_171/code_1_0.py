import datetime

def calculate_lc_refusal_deadline():
    """
    Calculates the deadline for a bank to send a refusal message for a Letter of Credit.
    """
    # --- Initial Parameters ---
    presentation_datetime = datetime.datetime(2020, 3, 30, 17, 1)
    bank_close_time = datetime.time(17, 0)
    # In Python's weekday(), Monday is 0 and Sunday is 6.
    non_banking_days = [5, 6]  # Saturday and Sunday

    print("Problem: A bank received documents at 05:01 pm on 30 March 2020. When is the latest they can send a refusal?")
    print("Rule: The bank has five banking days following the day of presentation (UCP 600, Article 16d).\n")

    # --- Step 1: Determine the effective day of presentation ---
    effective_presentation_date = presentation_datetime.date()
    print(f"1. Documents were presented on {presentation_datetime.strftime('%A, %d %B %Y')} at {presentation_datetime.strftime('%I:%M %p')}.")
    print(f"2. The bank's closing time is {bank_close_time.strftime('%I:%M %p')}.")

    if presentation_datetime.time() > bank_close_time:
        print("3. Since this is after closing time, the presentation is considered received on the next banking day.")
        next_day = effective_presentation_date + datetime.timedelta(days=1)
        while next_day.weekday() in non_banking_days:
            next_day += datetime.timedelta(days=1)
        effective_presentation_date = next_day
    
    print(f"4. The effective day of presentation is: {effective_presentation_date.strftime('%A, %d %B %Y')}.")

    # --- Step 2: Calculate the 5-day period ---
    examination_start_date = effective_presentation_date + datetime.timedelta(days=1)
    print(f"5. The 5-day examination period starts on the following day: {examination_start_date.strftime('%A, %d %B %Y')}.\n")
    
    print("6. Counting the five banking days (skipping non-banking days):")
    
    banking_days_counted = 0
    current_date = examination_start_date
    deadline_date = None

    while banking_days_counted < 5:
        if current_date.weekday() not in non_banking_days:
            banking_days_counted += 1
            print(f"   - Banking Day {banking_days_counted}: {current_date.strftime('%A, %d %B %Y')}")
        else:
            # This day is a weekend, so we skip it.
            pass
        
        if banking_days_counted == 5:
            deadline_date = current_date
            break
        
        current_date += datetime.timedelta(days=1)

    print(f"\nConclusion: The 5th banking day is {deadline_date.strftime('%d %B %Y')}.")
    print(f"Therefore, the bank must send the refusal message latest on {deadline_date.strftime('%d %B %Y')}.")

calculate_lc_refusal_deadline()
<<<D>>>