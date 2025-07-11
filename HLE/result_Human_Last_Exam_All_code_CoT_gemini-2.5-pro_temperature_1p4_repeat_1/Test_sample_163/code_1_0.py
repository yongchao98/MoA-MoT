def explain_surveillance_choice():
    """
    This function explains the rationale for selecting the best post-procedure
    surveillance program for a patient after SFA stenting.
    """

    # Clinical Case Summary
    patient_scenario = {
        "Procedure": "Angioplasty and stenting of the Superficial Femoral Artery (SFA)",
        "Primary Concern": "In-stent restenosis (re-narrowing)",
        "Highest Risk Period": "3 to 12 months post-procedure"
    }

    # Evaluation of Surveillance Methods
    surveillance_methods = {
        "Ankle-Brachial Index (ABI)": "Less sensitive for detecting restenosis. A significant narrowing must occur before ABI drops.",
        "Arterial Duplex Ultrasound": "The gold standard. Directly visualizes the stent and measures blood flow. Can detect restenosis early, before symptoms appear."
    }

    # Analysis of Answer Choices
    analysis = {
        "A": "Incorrect. Relies on ABI, which is not sensitive enough for stent surveillance.",
        "B": "Incorrect. Uses duplex but misses the critical 6-month follow-up.",
        "C": "Incorrect. Relies on ABI.",
        "D": "Correct. Uses the gold-standard duplex ultrasound at the most critical intervals (3, 6, and 12 months) when restenosis is most likely to develop.",
        "E": "Incorrect. Annual follow-up is too infrequent for the first year, which is the highest-risk period."
    }

    print("--- Rationale for Post-Procedure Surveillance ---")
    print(f"\n1. The patient underwent SFA stenting. The main goal of surveillance is to detect 'in-stent restenosis' early.")
    print(f"\n2. The best tool for this is Arterial Duplex Ultrasound. It is far more sensitive than an ABI for finding narrowing inside a stent.")
    print(f"   - An ABI only shows a problem after the narrowing is already severe.")
    print(f"   - A Duplex scan can find the problem much earlier.")

    print(f"\n3. Restenosis typically develops between 3 and 12 months after the procedure. Therefore, surveillance must be frequent during this first year.")

    print("\n--- Evaluating the Options ---")
    print(f"Choice A and C are not ideal because they only use ABI.")
    print(f"Choice E is not ideal because annual surveillance is too infrequent in the first year.")
    print(f"Choice B is better as it uses Duplex, but it omits the crucial 6-month check-up.")
    print(f"\nChoice D is the most appropriate program:")
    print(f"  - It uses the correct test: Arterial Duplex.")
    print(f"  - It schedules the test at the most important times to catch restenosis: 3 months, 6 months, 12 months, and then continues monitoring.")

explain_surveillance_choice()
print("\n<<<D>>>")