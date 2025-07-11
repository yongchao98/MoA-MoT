import datetime

def solve_lc_deadline():
    """
    Calculates the deadline for a bank to refuse a Letter of Credit presentation.
    """
    # Step 1: Define initial conditions
    reception_datetime = datetime.datetime(2020, 3, 30, 17, 1)  # 30 March 2020, 05:01 pm
    bank_close_time = datetime.time(17, 0)  # 5:00 pm
    days_to_examine = 5

    print(f"Documents were received on: {reception_datetime.strftime('%A, %d %B %Y at %I:%M %p')}")
    print(f"The bank's closing time is: {bank_close_time.strftime('%I:%M %p')}")

    # Step 2: Determine the effective day of presentation
    effective_presentation_date = reception_datetime.date()

    if reception_datetime.time() > bank_close_time:
        print("\nSince documents were received after closing, the presentation is considered received on the next banking day.")
        # Move to the next calendar day first
        effective_presentation_date += datetime.timedelta(days=1)
    
    # Ensure the effective date is a banking day (Mon-Fri). Monday is 0, Sunday is 6.
    while effective_presentation_date.weekday() >= 5: # 5 = Saturday, 6 = Sunday
        effective_presentation_date += datetime.timedelta(days=1)

    print(f"The effective day of presentation is: {effective_presentation_date.strftime('%A, %d %B %Y')}")
    print("-" * 50)

    # Step 3: Calculate the deadline
    # The counting period starts on the day FOLLOWING the presentation.
    print(f"The {days_to_examine}-day examination period begins on the day after presentation.")
    
    banking_days_counted = 0
    current_date = effective_presentation_date
    
    while banking_days_counted < days_to_examine:
        current_date += datetime.timedelta(days=1)
        # Check if the day is a banking day (Monday to Friday)
        if current_date.weekday() < 5:
            banking_days_counted += 1
            print(f"Banking Day {banking_days_counted}: {current_date.strftime('%A, %d %B %Y')}")
        else:
            print(f"Skipping non-banking day: {current_date.strftime('%A, %d %B %Y')}")
            
    deadline_date = current_date
    
    print("-" * 50)
    print("Conclusion:")
    print(f"The bank has {days_to_examine} banking days following the day of presentation ({effective_presentation_date.strftime('%d %B %Y')}) to send a refusal.")
    print(f"The final day to send the refusal message is the 5th banking day.")
    print(f"\nFinal Answer: The latest date is {deadline_date.strftime('%d %B %Y')}.")

if __name__ == "__main__":
    solve_lc_deadline()