def get_ampule_bud():
    """
    Determines and explains the Beyond-Use Date (BUD) for a single-dose ampule
    based on sterile compounding standards (USP <797>).
    """
    # --- Parameters of the Scenario ---
    container_type = "single dose container ampule"
    environment = "sterile environment"
    action = "puncture or opening"

    # --- Established Guideline ---
    # According to USP <797>, an opened ampule cannot be stored.
    bud = "immediately"
    disposal_rule = "Any unused portion must be discarded."

    # --- Reasoning ---
    reason = (
        "An ampule is a sealed glass container that is broken open and cannot be resealed. "
        "Once opened, it is an open system that is exposed to the environment. "
        "This exposure creates a risk of microbial contamination that cannot be mitigated, "
        "so the contents must be used without delay."
    )

    # --- Final Output ---
    print(f"The BUD for a {container_type} from the time of {action} in a {environment} is determined as follows:")
    print("\n--- Final Equation ---")
    print(f"BUD = Use {bud}")
    print("\n--- Explanation ---")
    print(f"Rule: {disposal_rule}")
    print(f"Reason: {reason}")

if __name__ == "__main__":
    get_ampule_bud()