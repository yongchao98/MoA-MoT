def calculate_weight_and_explain():
    """
    Calculates the newborn's weight in kilograms to confirm macrosomia
    and prints the supporting calculation.
    """
    pounds = 12
    ounces = 1
    
    # Conversion factors
    lb_to_kg = 0.453592
    oz_to_kg = 0.0283495
    
    # Calculate total weight in kg
    total_kg = (pounds * lb_to_kg) + (ounces * oz_to_kg)
    
    # Macrosomia is typically defined as a birth weight > 4.0 kg
    is_macrosomic = total_kg > 4.0
    
    print("This newborn has a weight of 12 lb 1 oz.")
    print("To confirm macrosomia, we convert this to kilograms:")
    print(f"({pounds} lb * {lb_to_kg} kg/lb) + ({ounces} oz * {oz_to_kg} kg/oz) = {total_kg:.2f} kg")
    print(f"A weight of {total_kg:.2f} kg is significantly above the 4.0 kg threshold for macrosomia.")
    print("\nThis macrosomia, combined with the respiratory distress and left lung density, points to a congenital diaphragmatic hernia, which is caused by a pleuroperitoneal membrane defect.")

calculate_weight_and_explain()