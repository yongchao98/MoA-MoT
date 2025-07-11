def get_bud_for_ampule():
    """
    Explains the Beyond-Use Date (BUD) for a single-dose ampule
    based on sterile compounding guidelines.
    """
    container_type = "single-dose container ampule"
    environment = "sterile environment"
    bud_duration = "Immediate Use"

    print(f"Container Type: {container_type}")
    print(f"Environment: Opened in a {environment}")
    print("-" * 40)
    print("Guideline Explanation:")
    print("An ampule is a glass container that cannot be resealed after it is broken open.")
    print("According to USP (United States Pharmacopeia) sterile compounding standards, once an ampule is opened, its contents must be used immediately.")
    print("Any unused portion must be discarded. It cannot be stored for any length of time.")
    print("-" * 40)
    print(f"The official Beyond-Use Date (BUD) is: {bud_duration}")

get_bud_for_ampule()