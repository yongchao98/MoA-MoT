def get_ampule_bud():
    """
    This function explains the Beyond-Use Date (BUD) for a single-dose ampule
    based on United States Pharmacopeia (USP) <797> guidelines.
    """

    container_type = "Single-dose ampule"
    environment = "Sterile environment (e.g., ISO Class 5)"

    print(f"Question: What is the BUD for a {container_type} opened in a {environment}?")
    print("-" * 70)
    print("Analysis based on USP <797> guidelines:")
    print("1. An ampule is a glass container that is broken open and cannot be resealed.")
    print("2. Once opened, it becomes an open-system container, highly susceptible to contamination.")
    print("3. Because of this high risk, the contents must be used without any delay.")
    print("4. Any unused portion must be discarded immediately and cannot be stored.")
    print("-" * 70)
    print("Final Answer:")
    print("The Beyond-Use Date (BUD) is: Immediate Use")

if __name__ == "__main__":
    get_ampule_bud()