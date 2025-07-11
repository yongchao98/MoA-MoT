import datetime

def solve_lc_deadline():
    """
    Calculates the deadline for a bank to send a refusal message for
    documents presented under a Letter of Credit (LC).
    """
    # --- Initial Data ---
    receipt_datetime = datetime.datetime(2020, 3, 30, 17, 1)
    bank_close_time = datetime.time(17, 0)
    # Per UCP 600, a bank has 5 banking days following the day of presentation
    days_to_examine = 5

    # --- Step 1: Determine the effective date of presentation ---
    print("Step 1: Determine the effective date of presentation.")
    print(f"Documents were received on: {receipt_datetime.strftime('%d %B %Y, at %I:%M %p')}")
    print(f"The bank's working hours end at: {bank_close_time.strftime('%I:%M %p')}")

    presentation_date = receipt_datetime.date()
    # If received after closing time, the presentation is deemed made on the next banking day.
    if receipt_datetime.time() > bank_close_time:
        print("Since the documents were received after banking hours, the presentation is considered to have been made on the next banking day.")
        # Find the next banking day
        next_day = presentation_date + datetime.timedelta(days=1)
        while next_day.weekday() >= 5:  # 5 is Saturday, 6 is Sunday
            next_day += datetime.timedelta(days=1)
        presentation_date = next_day
    
    print(f"Therefore, the effective date of presentation is: {presentation_date.strftime('%d %B %Y (%A)')}\n")

    # --- Step 2: Calculate the 5 banking days following presentation ---
    print("Step 2: Calculate the deadline (5 banking days following the day of presentation).")
    print(f"The calculation starts on the day after presentation, which is {(presentation_date + datetime.timedelta(days=1)).strftime('%d %B %Y')}.")
    
    current_date = presentation_date
    banking_days_counted = 0
    
    while banking_days_counted < days_to_examine:
        current_date += datetime.timedelta(days=1)
        # Check if the current day is a banking day (Monday to Friday)
        if current_date.weekday() < 5:
            banking_days_counted += 1
            print(f"Banking Day {banking_days_counted}: {current_date.strftime('%d %B %Y (%A)')}")
        else:
            print(f"Skipping {current_date.strftime('%d %B %Y (%A)')} as it is a non-banking day.")
            
    deadline_date = current_date
    print("\nStep 3: Conclusion.")
    print("According to UCP 600, the bank has a maximum of five banking days following the day of presentation to send a refusal message.")
    print(f"The fifth banking day is {deadline_date.strftime('%d %B %Y')}.")
    print("\nFinal Answer:")
    print(f"The bank should send the refusal message latest on {deadline_date.strftime('%d %B %Y')}.")

solve_lc_deadline()
<<<D>>>