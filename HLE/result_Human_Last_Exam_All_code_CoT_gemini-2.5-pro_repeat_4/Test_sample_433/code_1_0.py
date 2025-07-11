def get_bud_for_single_dose_ampule():
    """
    This function provides the Beyond-Use Date (BUD) for a single-dose ampule
    based on USP <797> guidelines for sterile compounding.
    """
    # The BUD for a single-dose container opened in an ISO Class 5 or cleaner environment.
    bud_hours = 12

    # Print the explanation and the result.
    print("The question is: What is the BUD for a single dose container ampule from the time of puncture in a sterile environment?")
    print("\nAccording to United States Pharmacopeia (USP) Chapter <797> guidelines:")
    print("A single-dose container (like an ampule), once opened or punctured within a sterile environment (ISO Class 5 or better), must be used within the following timeframe:")
    print(f"\nFinal Answer: {bud_hours} hours")
    print("\nNote: The opened ampule cannot be stored. Any unused portion must be discarded immediately after the initial use or at the end of the 12-hour period, whichever comes first.")

get_bud_for_single_dose_ampule()