def get_bud_for_ampule():
    """
    Calculates and explains the Beyond-Use Date (BUD) for a single-dose ampule.
    """
    container_type = "single-dose container ampule"
    environment = "sterile environment (ISO Class 5 or cleaner)"
    
    # According to USP Chapter <797> guidelines for sterile compounding:
    bud_in_sterile_env_hours = 12
    bud_in_suboptimal_env_hours = 1

    print(f"Container Type: {container_type}")
    print(f"Puncture Environment: {environment}")
    print("-" * 30)
    print("Based on USP <797> guidelines:")
    print(f"The Beyond-Use Date (BUD) for a {container_type} punctured in a {environment} is up to {bud_in_sterile_env_hours} hours.")
    print("\nNote: If opened in an environment of lower air quality (worse than ISO Class 5), the BUD is 1 hour.")
    print("Best practice for an open ampule is immediate use, with any unused portion discarded.")
    print("\nFinal Answer:")
    print(f"The maximum BUD is {bud_in_sterile_env_hours} hours.")

get_bud_for_ampule()