def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their magnetic field properties.

    The properties are:
    (A) driving a toroidal current
    (B) elongating flux surfaces / creating rotational transform
    (C) making the magnetic axis non-planar

    The classification logic is as follows:
    - Devices using only external coils with a planar axis are "Only B".
    - Devices using external coils with a non-planar axis are "Both B and C".
    - Devices driven by a plasma current are placed in the "A, B, and C" category,
      as this is the only category that includes property A.
    """
    
    classification = {
        "Only B": ["LHD"],
        "Both B and C": ["Wendelstein 7-X", "NCSX"],
        "A, B, and C": ["Tokamaks", "Reversed Field Pinches"]
    }

    print("Classification of Fusion Experiments:\n")
    
    # Print the results for each category
    category_b = classification["Only B"]
    print(f"Use Only Property B: {', '.join(category_b)}")
    
    category_bc = classification["Both B and C"]
    print(f"Use Both Properties B and C: {', '.join(category_bc)}")

    category_abc = classification["A, B, and C"]
    print(f"Use Properties A, B, and C: {', '.join(category_abc)}")

# Execute the function to print the classification
classify_fusion_experiments()