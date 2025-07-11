import datetime

def solve_lc_deadline():
    """
    Calculates the deadline for a bank to refuse a Letter of Credit presentation.
    """
    presentation_date = datetime.date(2020, 3, 30)
    banking_days_limit = 5
    
    print(f"The documents were presented on: {presentation_date.strftime('%A, %d %B %Y')}")
    print(f"According to UCP 600, the bank has a maximum of {banking_days_limit} banking days following the day of presentation to send a refusal.")
    print("The counting starts on the day after presentation. Saturdays and Sundays are non-banking days.\n")
    print("Let's count the banking days:")

    current_date = presentation_date
    banking_days_counted = 0
    
    while banking_days_counted < banking_days_limit:
        current_date += datetime.timedelta(days=1)
        # In Python's weekday(), Monday is 0 and Sunday is 6.
        # We consider days from Monday (0) to Friday (4) as banking days.
        if current_date.weekday() < 5:
            banking_days_counted += 1
            print(f"Banking Day {banking_days_counted}: {current_date.strftime('%A, %d %B %Y')}")
        else:
            print(f"Skipped       : {current_date.strftime('%A, %d %B %Y')} (Non-Banking Day)")

    deadline_date = current_date
    print(f"\nThe deadline is the end of the 5th banking day.")
    print(f"Final Answer: The latest date to send the refusal message is {deadline_date.strftime('%d %B %Y')}.")

solve_lc_deadline()