def classify_fusion_experiments():
    """
    Classifies fusion experiments based on the mechanisms they use to twist the magnetic field.

    The properties are:
    (A) Driving a toroidal current.
    (B) Elongating flux surfaces and making them rotate poloidally (3D shaping).
    (C) Making the magnetic axis non-planar.
    """

    classification = {
        "Uses only property B": [],
        "Uses both B and C": ["LHD", "Wendelstein 7-X"],
        "Uses A, B, and C": ["Tokamaks", "NCSX", "Reversed Field Pinches"]
    }

    print("Classification of Fusion Experiments:\n")

    category_b = "Uses only property B"
    experiments_b = classification[category_b]
    print(f"- {category_b}:")
    if experiments_b:
        for exp in experiments_b:
            print(f"  - {exp}")
    else:
        print("  - None")

    print("-" * 35)

    category_bc = "Uses both B and C"
    experiments_bc = classification[category_bc]
    print(f"- {category_bc}:")
    for exp in experiments_bc:
        print(f"  - {exp}")

    print("-" * 35)

    category_abc = "Uses A, B, and C"
    experiments_abc = classification[category_abc]
    print(f"- {category_abc}:")
    for exp in experiments_abc:
        print(f"  - {exp}")

classify_fusion_experiments()
<<<
Uses only property B:
  - None
-----------------------------------
Uses both B and C:
  - LHD
  - Wendelstein 7-X
-----------------------------------
Uses A, B, and C:
  - Tokamaks
  - NCSX
  - Reversed Field Pinches
>>>