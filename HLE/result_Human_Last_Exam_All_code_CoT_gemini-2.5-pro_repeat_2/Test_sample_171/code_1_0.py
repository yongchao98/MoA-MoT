import datetime

def solve():
    """
    Calculates the deadline for a bank to send a refusal message for a Letter of Credit presentation.
    """
    presentation_datetime = datetime.datetime(2020, 3, 30, 17, 1)
    bank_close_time = datetime.time(17, 0)
    
    print(f"Initial Presentation: {presentation_datetime.strftime('%d %B %Y, at %I:%M %p')}")
    print(f"Bank's Closing Hour: {bank_close_time.strftime('%I:%M %p')}\n")

    # Determine the official date of presentation
    presentation_date = presentation_datetime.date()
    if presentation_datetime.time() > bank_close_time:
        print("Presentation was made after banking hours.")
        # Find the next banking day
        effective_date = presentation_date + datetime.timedelta(days=1)
        while effective_date.weekday() >= 5: # 5 is Saturday, 6 is Sunday
            effective_date += datetime.timedelta(days=1)
        print(f"The presentation is considered received on the next banking day: {effective_date.strftime('%A, %d %B %Y')}\n")
    else:
        effective_date = presentation_date
        print(f"The presentation is considered received on: {effective_date.strftime('%A, %d %B %Y')}\n")
        
    # The examination period starts on the day FOLLOWING the presentation day.
    # We need to find 5 banking days after the effective_date.
    
    print("Calculating the 5 banking days for examination:")
    
    current_date = effective_date
    banking_days_to_count = 5
    days_counted = 0
    
    while days_counted < banking_days_to_count:
        current_date += datetime.timedelta(days=1)
        # Check if the day is a banking day (Monday to Friday)
        if current_date.weekday() < 5:
            days_counted += 1
            day_type = "Banking Day"
            print(f"Day {days_counted}: {current_date.strftime('%A, %d %B %Y')} ({day_type})")
        else:
            day_type = "Non-Banking Day"
            print(f"       {current_date.strftime('%A, %d %B %Y')} ({day_type} - Skipped)")

    deadline_date = current_date
    print(f"\nThe bank has until the end of the 5th banking day to send a refusal message.")
    print(f"Final Deadline: Latest on {deadline_date.strftime('%d %B %Y')}")

solve()
<<<D>>>