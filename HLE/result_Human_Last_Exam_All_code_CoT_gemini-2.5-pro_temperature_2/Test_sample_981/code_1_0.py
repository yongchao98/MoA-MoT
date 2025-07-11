def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their magnetic field twisting methods
    and prints the results.
    """

    # A: driving a toroidal current
    # B: elongating the flux surfaces and making them rotate poloidally
    # C: making the magnetic axis non-planar (or having other 3D field effects like ripple)

    classification = {
        "using Only property B": ["LHD", "Wendelstein 7-X"],
        "using both properties B and C": ["NCSX"],
        "using properties A, B, and C": ["Tokamaks", "Reversed Field Pinches"]
    }

    print("Classification of Fusion Experiments by Magnetic Field Twisting Method:")
    print("----------------------------------------------------------------------")
    print("Property Definitions:")
    print("  (A) Driving a toroidal current")
    print("  (B) Elongating flux surfaces and making them rotate poloidally")
    print("  (C) Making the magnetic axis non-planar\n")

    for category, experiments in classification.items():
        # Format the output line with the category name
        line = f"Experiments {category}:"
        print(line)
        # Print each experiment in that category
        for experiment in experiments:
            print(f"  - {experiment}")
        print() # Add a blank line for readability

if __name__ == '__main__':
    classify_fusion_experiments()