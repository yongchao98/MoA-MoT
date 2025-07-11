import datetime

def calculate_lc_deadline():
    """
    Calculates the deadline for a bank to send a refusal message for a Letter of Credit.
    """
    # Step 1: Define initial parameters
    presentation_datetime = datetime.datetime(2020, 3, 30, 17, 1)
    bank_close_time = datetime.time(17, 0)
    non_banking_days = [5, 6]  # Saturday = 5, Sunday = 6
    days_to_examine = 5

    print(f"Documents received on: {presentation_datetime.strftime('%A, %d %B %Y at %I:%M %p')}")
    print(f"Bank's closing time: {bank_close_time.strftime('%I:%M %p')}\n")

    # Step 2: Determine the effective day of presentation
    effective_presentation_date = presentation_datetime.date()
    if presentation_datetime.time() > bank_close_time:
        print("Since documents were received after closing time, the presentation is deemed made on the next banking day.")
        # Find the next banking day
        next_day = effective_presentation_date + datetime.timedelta(days=1)
        while next_day.weekday() in non_banking_days:
            next_day += datetime.timedelta(days=1)
        effective_presentation_date = next_day
    
    print(f"Effective date of presentation: {effective_presentation_date.strftime('%A, %d %B %Y')}\n")

    # Step 3: Calculate the 5 banking days following the effective presentation date
    print(f"The examination period is {days_to_examine} banking days starting from the day after presentation.")
    
    current_date = effective_presentation_date
    banking_days_counted = 0
    
    while banking_days_counted < days_to_examine:
        current_date += datetime.timedelta(days=1)
        # Check if it's a banking day
        if current_date.weekday() not in non_banking_days:
            banking_days_counted += 1
            print(f"Day {banking_days_counted}: {current_date.strftime('%A, %d %B %Y')}")
        else:
            print(f"Skip: {current_date.strftime('%A, %d %B %Y')} (Non-banking day)")
            
    deadline_date = current_date
    print(f"\nThe bank must send the refusal message at the latest on this day.")
    print(f"Final Deadline: {deadline_date.strftime('%d %B %Y')}")

calculate_lc_deadline()
<<<D>>>