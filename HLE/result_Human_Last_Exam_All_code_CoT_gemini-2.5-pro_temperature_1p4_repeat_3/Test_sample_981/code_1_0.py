def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their magnetic field twisting properties.

    The properties are:
    A: Driving a toroidal current.
    B: Elongating flux surfaces and making them rotate poloidally (rotational transform).
    C: Making the magnetic axis non-planar.

    The classification logic is as follows:
    - Devices with a non-planar axis (C) are stellarators.
      - If they also rely on a significant driven current (A), they are A, B, and C (e.g., hybrid stellarators).
      - If they do not rely on a driven current, they are B and C (e.g., pure stellarators).
    - Devices with a planar axis (not C) are tokamaks or RFPs.
      - They are classified as "Only B", interpreting the current (A) as the means to achieve the twist (B).
    """

    # Categories based on the properties A, B, C
    only_B = {
        "name": "Only property B",
        "experiments": ["Tokamaks", "Reversed Field Pinches"]
    }

    B_and_C = {
        "name": "Properties B and C",
        "experiments": ["LHD", "Wendelstein 7-X"]
    }

    A_B_and_C = {
        "name": "Properties A, B, and C",
        "experiments": ["NCSX"]
    }

    # Print the results in a structured format
    print("Classification of Fusion Experiments:\n")

    print(f"{only_B['name']}:")
    for exp in only_B['experiments']:
        print(f"- {exp}")

    print(f"\n{B_and_C['name']}:")
    for exp in B_and_C['experiments']:
        print(f"- {exp}")

    print(f"\n{A_B_and_C['name']}:")
    for exp in A_B_and_C['experiments']:
        print(f"- {exp}")

if __name__ == "__main__":
    classify_fusion_experiments()
    # The final answer is the derived classification.
    # The code prints the full classification as requested.
    # We will format the final summary answer as a string.
    final_answer = "Only B: Tokamaks, Reversed Field Pinches; B and C: LHD, Wendelstein 7-X; A, B, and C: NCSX"
    print(f"\n<<<{final_answer}>>>")
