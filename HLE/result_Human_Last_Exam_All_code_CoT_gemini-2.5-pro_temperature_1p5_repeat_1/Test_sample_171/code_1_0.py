import datetime

def calculate_lc_refusal_deadline():
    """
    Calculates the deadline for a bank to refuse documents under a Letter of Credit
    based on UCP 600 rules.
    """
    # --- Given Information ---
    presentation_datetime = datetime.datetime(2020, 3, 30, 17, 1)  # 30 March 2020, 05:01 pm
    bank_close_time = datetime.time(17, 0)  # 5:00 pm
    days_for_examination = 5
    # Note: In Python's weekday(), Monday is 0 and Sunday is 6.
    # Non-banking days are Saturday (5) and Sunday (6).

    print("--- Step 1: Determine the Official Day of Presentation ---")
    effective_date = presentation_datetime.date()
    print(f"Initial presentation date and time: {presentation_datetime.strftime('%A, %d %B %Y, %I:%M %p')}")
    print(f"Bank's closing time: {bank_close_time.strftime('%I:%M %p')}")

    # If presentation is after closing time, it's considered made on the next day.
    if presentation_datetime.time() > bank_close_time:
        print("Presentation was after closing time. It is deemed received on the next business day.")
        effective_date += datetime.timedelta(days=1)
    
    # Ensure the effective date is a banking day (not Sat/Sun)
    while effective_date.weekday() >= 5:
        print(f"Adjusting: {effective_date.strftime('%A, %d %B %Y')} is a non-banking day. Moving to the next day.")
        effective_date += datetime.timedelta(days=1)

    print(f"\nThe official 'Day of Presentation' is: {effective_date.strftime('%A, %d %B %Y')}")
    print("\n" + "-"*50)

    print("--- Step 2: Calculate the 5 Banking Day Deadline ---")
    print("The examination period starts the day AFTER the Day of Presentation.")
    
    current_date = effective_date
    banking_days_counted = 0
    
    while banking_days_counted < days_for_examination:
        current_date += datetime.timedelta(days=1)
        # Check if the day is a banking day (Monday to Friday)
        if current_date.weekday() < 5:
            banking_days_counted += 1
            print(f"Banking Day {banking_days_counted}: {current_date.strftime('%A, %d %B %Y')}")
        else:
            print(f"Skipped       : {current_date.strftime('%A, %d %B %Y')} (Non-banking day)")

    deadline_date = current_date
    print("\n" + "-"*50)

    print("--- Step 3: Final Deadline ---")
    print(f"The calculation starts from {effective_date.strftime('%d %B %Y')} and counts {days_for_examination} banking days.")
    print("The final day is the deadline for sending the refusal message.")
    print(f"\nFinal Answer: The latest date is {deadline_date.strftime('%d %B %Y')}.")

if __name__ == '__main__':
    calculate_lc_refusal_deadline()