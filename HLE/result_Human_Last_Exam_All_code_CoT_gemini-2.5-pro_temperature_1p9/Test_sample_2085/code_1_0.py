def find_french_equivalent():
    """
    This function demonstrates the linguistic transformation of the name
    "Odobeccus" into its French equivalent, "Odobec".
    """
    original_name = "Odobeccus"
    step1_name = "Odobecc"
    final_name = "Odobec"

    print(f"Finding the French equivalent for '{original_name}':")
    print(f"Step 1: Drop the Latin nominative suffix '-us'.")
    print(f"{original_name} -> {step1_name}")
    print("\nStep 2: Simplify the final double consonant 'cc' to 'c', a common simplification in French.")
    print(f"{step1_name} -> {final_name}")
    print(f"\nThe resulting French equivalent is: {final_name}")

find_french_equivalent()