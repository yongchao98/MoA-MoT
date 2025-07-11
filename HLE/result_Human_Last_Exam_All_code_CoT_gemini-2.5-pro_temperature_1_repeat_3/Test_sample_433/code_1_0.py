import sys

def get_bud_for_single_dose_container(container_type, environment):
    """
    Determines and prints the Beyond-Use Date (BUD) for a single-dose container
    based on USP <797> guidelines.

    Args:
        container_type (str): The type of container, e.g., "ampule" or "vial".
        environment (str): The environment where the container is opened, e.g., "sterile".
    """
    print(f"Calculating BUD for a single-dose '{container_type}' opened in a '{environment}' environment.")
    print("="*70)

    if container_type.lower() == "ampule":
        # An ampule is a glass container that is broken open and cannot be resealed.
        # It is considered an open system once the top is snapped off.
        storage_time_hours = 0
        recommendation = "Immediate Use Only"
        
        print(f"Rule: An opened ampule cannot be stored.")
        print(f"Explanation: Once an ampule is opened, it is exposed to the environment and cannot be resealed.")
        print(f"Calculated Storage Time: {storage_time_hours} hours.")
        print(f"Final BUD Guideline: {recommendation}. The contents must be used immediately, and any unused portion must be discarded.")
        
    elif container_type.lower() == "vial":
        # A single-dose vial has a rubber stopper that can be punctured.
        if environment.lower() == "sterile":
            storage_time_hours = 12
            print(f"Rule: A single-dose vial punctured in an ISO Class 5 or cleaner environment.")
            print(f"Final BUD Guideline: May be used for up to {storage_time_hours} hours after puncture.")
        else:
            storage_time_hours = 1
            print(f"Rule: A single-dose vial punctured in an environment worse than ISO Class 5.")
            print(f"Final BUD Guideline: Must be used within {storage_time_hours} hour or discarded.")
    else:
        print(f"Error: Unknown container type '{container_type}'.", file=sys.stderr)

# User's specific case
get_bud_for_single_dose_container(container_type="ampule", environment="sterile")
