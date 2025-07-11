def get_ampule_bud_information():
    """
    Provides the Beyond-Use Date (BUD) information for a single-dose ampule.
    """
    rule = "Immediate use"
    explanation = "According to USP <797> guidelines, once a single-dose ampule is opened, it is exposed to the environment and cannot be resealed."
    disposal_instruction = "The contents must be used immediately, and any unused portion must be discarded."
    conclusion = "There is no permissible storage time, so the Beyond-Use Date is effectively immediate."

    print(f"BUD Information for an opened Single-Dose Ampule:\n")
    print(f"Guideline: {rule}")
    print(f"Reason: {explanation}")
    print(f"Action: {disposal_instruction}")
    print(f"Conclusion: {conclusion}")

get_ampule_bud_information()