import datetime

def solve_lc_deadline():
    """
    Calculates the deadline for a bank to refuse a Letter of Credit presentation.
    """
    # Initial presentation details
    presentation_datetime = datetime.datetime(2020, 3, 30, 17, 1)
    bank_close_time = datetime.time(17, 0)
    
    print(f"Initial presentation date and time: {presentation_datetime.strftime('%A, %d %B %Y, %I:%M %p')}")
    print(f"Bank's closing hour: {bank_close_time.strftime('%I:%M %p')}\n")

    # Step 1: Determine the official date of presentation
    presentation_date = presentation_datetime.date()
    effective_presentation_date = presentation_date

    print("Step 1: Determine the official presentation date.")
    if presentation_datetime.time() > bank_close_time:
        print("Presentation was received after banking hours.")
        # Find the next banking day
        day_increment = 1
        while True:
            next_day = presentation_date + datetime.timedelta(days=day_increment)
            # Monday is 0, Sunday is 6. We skip Saturday (5) and Sunday (6).
            if next_day.weekday() not in [5, 6]:
                effective_presentation_date = next_day
                break
            day_increment += 1
    
    print(f"The official day of presentation is considered to be: {effective_presentation_date.strftime('%A, %d %B %Y')}\n")

    # Step 2: Calculate the 5 banking days following the presentation date
    print("Step 2: Calculate the 5 banking days allowed for refusal, starting from the day AFTER presentation.")
    
    banking_days_to_count = 5
    days_counted = 0
    current_date = effective_presentation_date
    final_date = None

    start_date_str = effective_presentation_date.strftime('%d %B %Y')

    while days_counted < banking_days_to_count:
        current_date += datetime.timedelta(days=1)
        if current_date.weekday() not in [5, 6]:
            days_counted += 1
            print(f"Banking Day {days_counted}: {current_date.strftime('%A, %d %B %Y')}")
            if days_counted == banking_days_to_count:
                final_date = current_date

    final_date_str = final_date.strftime('%d %B %Y')

    print("\n--- Conclusion ---")
    print("The bank has until the end of the fifth banking day to send a refusal message.")
    print(f"The final equation is: {start_date_str} + {banking_days_to_count} banking days = {final_date_str}")
    print(f"\nTherefore, the latest date for refusal is {final_date.strftime('%d %B %Y')}.")

solve_lc_deadline()
<<<D>>>