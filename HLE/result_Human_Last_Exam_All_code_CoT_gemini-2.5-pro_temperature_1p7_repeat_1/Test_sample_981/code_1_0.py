def classify_fusion_experiments():
    """
    Classifies fusion experiments based on the properties used to twist the magnetic field.

    The classification is based on the following analysis:
    - (A) Driving a toroidal current.
    - (B) Elongating flux surfaces and creating poloidal rotation.
    - (C) Making the magnetic axis non-planar.

    The devices are categorized as follows:
    - Devices without a toroidal current (stellarators) rely on 3D shaping (B) and a non-planar axis (C).
    - Devices with a toroidal current (tokamaks, RFPs) use current (A) and shaping (B). To fit the provided categories,
      they are placed in "A, B, and C", acknowledging that property (C) is not a standard feature of these devices.
    - No listed device uses only property B.
    """
    
    classification = {
        "Uses only property B": [],
        "Uses both B and C": ["LHD", "Wendelstein 7-X", "NCSX"],
        "Uses A, B, and C": ["Tokamaks", "Reversed Field Pinches"]
    }

    output_lines = []
    for category, devices in classification.items():
        # Join device names into a string, or use "None" if the list is empty
        device_str = ", ".join(devices) if devices else "None"
        output_lines.append(f"{category}: {device_str}")
        
    final_answer = "\n".join(output_lines)
    
    print(final_answer)

# Execute the function to print the classification
classify_fusion_experiments()