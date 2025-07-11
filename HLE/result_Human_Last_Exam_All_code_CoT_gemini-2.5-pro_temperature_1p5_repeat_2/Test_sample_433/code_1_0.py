def get_bud_information():
    """
    Provides Beyond-Use Date (BUD) information for single-dose containers
    based on USP <797> guidelines.
    """
    bud_rules = {
        "ampule": {
            "name": "Single-Dose Ampule",
            "environment": "Any (including sterile ISO Class 5)",
            "bud": "Immediate use",
            "explanation": "Once an ampule is opened, it cannot be resealed. It is an open system exposed to the environment, so its contents must be used immediately, and any unused portion must be discarded."
        },
        "single_dose_vial_iso5": {
            "name": "Single-Dose Vial",
            "environment": "Sterile (ISO Class 5 or cleaner)",
            "bud_hours": 12,
            "explanation": "A single-dose vial punctured in a sterile environment can be used for up to 12 hours, as its rubber stopper can be disinfected and reseals after puncture."
        }
    }

    # Answer the specific question about the ampule
    ampule_info = bud_rules["ampule"]
    print("--- BUD for a Single-Dose Ampule ---")
    print(f"Container Type: {ampule_info['name']}")
    print(f"Punctured in Environment: {ampule_info['environment']}")
    print(f"Beyond-Use Date (BUD): {ampule_info['bud']}")
    print(f"Reason: {ampule_info['explanation']}")
    print("-" * 35)

    # Provide context by showing the rule for a single-dose vial
    vial_info = bud_rules["single_dose_vial_iso5"]
    print("\n--- For Comparison: BUD for a Single-Dose Vial ---")
    print(f"Container Type: {vial_info['name']}")
    print(f"Punctured in Environment: {vial_info['environment']}")
    # Output the number from the rule as requested
    bud_duration_hours = vial_info['bud_hours']
    print(f"Beyond-Use Date (BUD): Up to {bud_duration_hours} hours after initial puncture.")
    print(f"Reason: {vial_info['explanation']}")
    print("-" * 50)


if __name__ == "__main__":
    get_bud_information()