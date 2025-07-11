def find_treatment_recommendation():
    """
    This function prints key sentences from the provided text to answer the user's question.
    """

    first_line_therapy = "Current National Comprehensive Cancer Network (NCCN) guidelines with category 1 evidence recommend THP as the preferred first-line regimen for HER2+ MBC."
    second_line_therapy = "Upon disease progression, TDXd is recommended in the second-line setting."
    third_line_therapy = "Although the optimal sequence for third-line therapy is not known, NCCN guidelines recommend tucatinib, trastuzumab, and capecitabine as the preferred regimen."
    replacement_of_old_therapy = "In 2022, TDXd replaced T-DM1 as second-line therapy on the basis of superior results reported in the DESTINY-Breast03 trial."

    print("Relevant information from the text regarding treatment progression for HER2+ MBC:")
    print("--------------------------------------------------------------------------------")
    print("First-Line (The treatment the patient has already had):")
    print(f'"{first_line_therapy}" (THP includes a taxane and trastuzumab.)')
    print("\nSecond-Line (After progression on first-line):")
    print(f'"{second_line_therapy}" (TDXd is fam-trastuzumab deruxtecan.)')
    print("\nReasoning for Replacement:")
    print(f'"{replacement_of_old_therapy}" (This rules out T-DM1 / trastuzumab emtansine as the current recommended second-line option.)')
    print("\nThird-Line (A later option):")
    print(f'"{third_line_therapy}" (This rules out tucatinib, trastuzumab, and capecitabine as the second-line option.)')
    print("--------------------------------------------------------------------------------")
    print("Conclusion: The text indicates that Fam-trastuzumab deruxtecan (TDXd) is the recommended second-line therapy.")

find_treatment_recommendation()