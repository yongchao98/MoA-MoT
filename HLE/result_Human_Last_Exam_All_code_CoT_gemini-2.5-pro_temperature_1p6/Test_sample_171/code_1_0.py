import datetime

def calculate_lc_refusal_deadline():
    """
    Calculates the deadline for a bank to send a refusal message for a Letter of Credit.
    """
    receipt_datetime_str = "2020-03-30 17:01"
    receipt_datetime = datetime.datetime.strptime(receipt_datetime_str, "%Y-%m-%d %H:%M")
    bank_closing_time = datetime.time(17, 0)
    days_to_examine = 5

    # --- Step 1: Determine the effective date of presentation ---
    print("Step 1: Determine the effective date of presentation.")
    print(f"Documents received on: {receipt_datetime.strftime('%A, %d %B %Y at %I:%M %p')}")
    print(f"Bank closing time: {bank_closing_time.strftime('%I:%M %p')}")

    presentation_date = receipt_datetime.date()
    # If received after closing time, presentation is considered for the next banking day.
    if receipt_datetime.time() > bank_closing_time:
        print("This is after banking hours. Presentation is deemed made on the next banking day.")
        presentation_date += datetime.timedelta(days=1)
    
    # Ensure the presentation date is a banking day (Mon-Fri)
    while presentation_date.weekday() >= 5: # 5 = Saturday, 6 = Sunday
        presentation_date += datetime.timedelta(days=1)

    print(f"Effective presentation date: {presentation_date.strftime('%A, %d %B %Y')}\n")

    # --- Step 2: Apply UCP 600 and count the five banking days ---
    print(f"Step 2: Calculate the deadline based on UCP 600.")
    print(f"The rule allows a maximum of {days_to_examine} banking days following the day of presentation.")
    print("The counting starts on the day after the effective presentation date.\n")
    
    print("Counting the five banking days:")
    banking_days_counted = 0
    current_date = presentation_date
    while banking_days_counted < days_to_examine:
        current_date += datetime.timedelta(days=1)
        # Weekday check: Monday is 0, Sunday is 6. We only count days < 5.
        if current_date.weekday() < 5:
            banking_days_counted += 1
            print(f"   Banking Day {banking_days_counted}: {current_date.strftime('%A, %d %B %Y')}")
        else:
            print(f"   Skipping: {current_date.strftime('%A, %d %B %Y')} (weekend)")
            
    deadline_date = current_date
    
    # --- Step 3: Announce the final deadline ---
    print("\nStep 3: State the final deadline.")
    print(f"The fifth and final banking day is {deadline_date.strftime('%A, %d %B %Y')}.")
    print(f"Therefore, the bank must send the refusal message latest on {deadline_date.strftime('%d %B %Y')}.")

if __name__ == "__main__":
    calculate_lc_refusal_deadline()