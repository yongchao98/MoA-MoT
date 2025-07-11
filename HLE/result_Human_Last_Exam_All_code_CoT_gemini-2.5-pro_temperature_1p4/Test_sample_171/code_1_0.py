import datetime

def solve_lc_deadline():
    """
    Calculates the deadline for a bank to send a refusal message for a Letter of Credit.
    """
    # --- Initial Data ---
    presentation_datetime_str = "2020-03-30 17:01"
    bank_close_time_str = "17:00"
    
    presentation_dt = datetime.datetime.strptime(presentation_datetime_str, "%Y-%m-%d %H:%M")
    bank_close_time = datetime.datetime.strptime(bank_close_time_str, "%H:%M").time()
    
    # weekday(): Monday is 0 and Sunday is 6. Saturday is 5.
    non_banking_days_weekday = [5, 6] 

    print("Step 1: Determine the official 'Day of Presentation'.")
    print(f"Documents were physically received on: {presentation_dt.strftime('%A, %d %B %Y at %I:%M %p')}")
    print(f"The bank's working hours end at: {bank_close_time.strftime('%I:%M %p')}")

    # --- Determine Effective Presentation Date ---
    effective_presentation_date = presentation_dt.date()
    if presentation_dt.time() > bank_close_time:
        print("\nSince the documents were received after banking hours, the presentation is deemed to be on the next banking day.")
        # Move to the next day and find the first available banking day
        effective_presentation_date += datetime.timedelta(days=1)
        while effective_presentation_date.weekday() in non_banking_days_weekday:
            effective_presentation_date += datetime.timedelta(days=1)
    
    print(f"The official 'Day of Presentation' is: {effective_presentation_date.strftime('%A, %d %B %Y')}")

    print("\nStep 2: Calculate the 5 banking-day examination period.")
    print("The examination period starts on the day *after* the Day of Presentation.")

    # --- Calculate the Deadline ---
    current_date = effective_presentation_date
    banking_days_counted = 0
    days_list = []

    while banking_days_counted < 5:
        # Move to the next calendar day
        current_date += datetime.timedelta(days=1)
        # Check if it's a banking day
        if current_date.weekday() not in non_banking_days_weekday:
            banking_days_counted += 1
            days_list.append(current_date)
            
    deadline_date = days_list[-1]

    print("\nThe final calculation is:")
    print(f"Official Presentation Day ({effective_presentation_date.strftime('%d %b')}) + Banking Day 1 = {days_list[0].strftime('%d %b')}")
    print(f"+ Banking Day 2 = {days_list[1].strftime('%d %b')}")
    print(f"+ Banking Day 3 = {days_list[2].strftime('%d %b')}")
    print(f"+ Banking Day 4 = {days_list[3].strftime('%d %b')}")
    print(f"+ Banking Day 5 = {days_list[4].strftime('%d %b')}")

    print("\nStep 3: State the final deadline.")
    print(f"The refusal message must be sent at the latest by the close of business on the 5th banking day.")
    print(f"The deadline is: {deadline_date.strftime('%d %B %Y')}")


solve_lc_deadline()
<<<D>>>