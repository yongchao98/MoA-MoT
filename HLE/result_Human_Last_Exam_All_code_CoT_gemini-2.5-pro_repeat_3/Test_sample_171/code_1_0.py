import datetime

def calculate_refusal_deadline():
    """
    Calculates the deadline for a bank to send a refusal message for a Letter of Credit.
    """
    presentation_datetime = datetime.datetime(2020, 3, 30, 17, 1)
    bank_close_time = datetime.time(17, 0, 0)

    print(f"Step 1: Determine the official presentation date.")
    print(f"Documents presented on: {presentation_datetime.strftime('%A, %d %B %Y at %I:%M %p')}")
    print(f"Bank closes at: {bank_close_time.strftime('%I:%M %p')}")

    # If presentation is after banking hours, it's deemed to be on the next banking day.
    presentation_date = presentation_datetime.date()
    if presentation_datetime.time() > bank_close_time:
        next_day = presentation_date + datetime.timedelta(days=1)
        # Ensure the next day is a banking day
        while next_day.weekday() >= 5: # 5 is Saturday, 6 is Sunday
            next_day += datetime.timedelta(days=1)
        official_presentation_date = next_day
        print(f"Since presentation was after hours, it is considered received on the next banking day.")
    else:
        official_presentation_date = presentation_date

    print(f"Official presentation date: {official_presentation_date.strftime('%A, %d %B %Y')}\n")

    print("Step 2: Count five banking days following the presentation date.")
    # The counting period starts on the day after the official presentation date.
    current_date = official_presentation_date + datetime.timedelta(days=1)
    banking_days_to_count = 5
    banking_days_counted = 0
    
    while banking_days_counted < banking_days_to_count:
        # Check if the day is a weekday (Monday=0, Sunday=6)
        if current_date.weekday() < 5:
            banking_days_counted += 1
            print(f"Banking Day {banking_days_counted}: {current_date.strftime('%A, %d %B %Y')}")
        else:
            print(f"Skipping non-banking day: {current_date.strftime('%A, %d %B %Y')}")
        
        # If we have found the last day, store it and break
        if banking_days_counted == banking_days_to_count:
            deadline_date = current_date
            break
        
        current_date += datetime.timedelta(days=1)

    print(f"\nFinal Answer: The refusal message must be sent latest on the fifth banking day.")
    print(f"The deadline is latest on {deadline_date.strftime('%d %B %Y')}.")

calculate_refusal_deadline()
<<<D>>>