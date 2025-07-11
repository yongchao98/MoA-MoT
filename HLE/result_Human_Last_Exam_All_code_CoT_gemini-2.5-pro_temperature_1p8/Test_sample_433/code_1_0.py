import sys

def get_single_dose_bud(container_type, environment):
    """
    Determines the Beyond-Use Date (BUD) for single-dose containers based on USP <797> guidelines.

    Args:
        container_type (str): The type of container ('ampule' or 'vial').
        environment (str): The environment where the container is opened ('sterile' for ISO 5 or 'non-sterile').
    """
    container_type = container_type.lower()
    environment = environment.lower()

    print(f"Calculating BUD for a single-dose '{container_type}' opened in a '{environment}' environment...")
    print("---")

    if container_type == "ampule":
        bud = "Immediate Use"
        explanation = ("According to USP <797> guidelines, once a glass ampule is opened, "
                       "it cannot be resealed. Therefore, its contents must be used immediately. "
                       "Any unused portion must be discarded and cannot be stored for any period.")
        print(f"Rule: {explanation}")
        print("---")
        print(f"Final BUD: {bud}")
        # We use a special return format for the final answer expected by the system.
        return f"<<<{bud}>>>"

    elif container_type == "vial":
        if environment == "sterile":
            bud = "12 hours"
            explanation = ("A single-dose vial (with a rubber stopper) entered in an ISO Class 5 "
                           "or cleaner environment may be used for up to 12 hours after the initial puncture.")
        else:
            bud = "1 hour"
            explanation = ("A single-dose vial entered in an environment worse than ISO Class 5 "
                           "must be used within 1 hour or discarded.")
        print(f"Rule: {explanation}")
        print("---")
        print(f"Final BUD: {bud}")
        return f"<<<{bud}>>>"
    else:
        print(f"Error: Container type '{container_type}' not recognized. Use 'ampule' or 'vial'.")
        return "<<<Invalid Input>>>"


if __name__ == "__main__":
    # This specifically addresses the user's question.
    # The return value from the function is captured but not printed here,
    # as the final answer is meant to be at the very end of the response.
    get_single_dose_bud(container_type="ampule", environment="sterile")