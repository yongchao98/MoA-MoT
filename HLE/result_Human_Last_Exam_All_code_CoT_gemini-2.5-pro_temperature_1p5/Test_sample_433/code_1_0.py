import datetime

def calculate_bud_single_dose_ampule():
    """
    Calculates and explains the Beyond-Use Date (BUD) for a single-dose ampule
    punctured in a sterile environment, based on USP <797> guidelines.
    """
    # According to USP <797>, the BUD for a single-dose container (like an ampule)
    # punctured or opened in an ISO Class 5 or cleaner environment is 12 hours.
    bud_in_hours = 12

    # Get the current time to simulate the moment of puncture.
    puncture_time = datetime.datetime.now()

    # Calculate the exact expiration time.
    expiration_time = puncture_time + datetime.timedelta(hours=bud_in_hours)

    # --- Outputting the results ---
    print("Based on USP <797> guidelines for a sterile (ISO Class 5 or cleaner) environment:")
    print(f"The Beyond-Use Date (BUD) for a single-dose ampule is {bud_in_hours} hours after puncture.")
    print("\n--- Example Calculation ---")
    print(f"Time of Puncture : {puncture_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"BUD Duration     : + {bud_in_hours} hours")
    print(f"Must Be Used By  : {expiration_time.strftime('%Y-%m-%d %H:%M:%S')}")

if __name__ == "__main__":
    calculate_bud_single_dose_ampule()