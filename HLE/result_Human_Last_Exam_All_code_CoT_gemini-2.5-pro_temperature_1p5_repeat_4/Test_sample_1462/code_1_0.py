def analyze_newborn_weight(pounds, ounces):
    """
    Converts newborn weight from pounds and ounces to kilograms and analyzes it.
    """
    # Conversion factors
    ounces_per_pound = 16
    kg_per_pound = 0.453592

    # --- Calculation ---
    # Convert total weight to pounds
    total_pounds = pounds + ounces / ounces_per_pound

    # Convert pounds to kilograms
    total_kg = total_pounds * kg_per_pound

    # --- Outputting the results ---
    print("Analyzing the newborn's weight to understand the clinical context.")
    print(f"Weight: {pounds} lbs {ounces} oz.")
    
    # Show the equation steps
    print("\nCalculation steps to convert to kilograms:")
    print(f"1. Total pounds = {pounds} + ({ounces} / {ounces_per_pound}) = {total_pounds:.4f} lbs")
    print(f"2. Total kg = {total_pounds:.4f} lbs * {kg_per_pound} kg/lb = {total_kg:.2f} kg")

    print(f"\nThe newborn's weight is approximately {total_kg:.2f} kg.")
    print("A typical newborn weighs around 3.5 kg. A weight over 4.5 kg is considered severe macrosomia.")
    print("This finding is a major clue, often associated with maternal diabetes, which is a risk factor for certain congenital defects.")

# Given data from the problem
newborn_pounds = 12
newborn_ounces = 1

analyze_newborn_weight(newborn_pounds, newborn_ounces)