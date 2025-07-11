def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their method of twisting the magnetic field.

    The properties are:
    A: Driving a toroidal current.
    B: Elongating flux surfaces and making them rotate poloidally.
    C: Making the magnetic axis non-planar.
    """

    # The classification is stored in a dictionary.
    # The keys are the categories of properties used.
    # The values are lists of the experiments that fall into each category.
    classification = {
        "Only property B": [],
        "Both B and C": [
            "LHD",
            "Wendelstein 7-X",
            "NCSX"
        ],
        "A, B, and C": [
            "Tokamaks",
            "Reversed Field Pinches"
        ]
    }

    print("Classification of Fusion Experiments by Field Twisting Properties:")
    print("-" * 60)
    print("Property A: Driving a toroidal current")
    print("Property B: Elongating flux surfaces and poloidal rotation")
    print("Property C: Non-planar magnetic axis")
    print("-" * 60 + "\n")


    for category, experiments in classification.items():
        print(f"Experiments that use {category}:")
        if experiments:
            for exp in experiments:
                print(f"- {exp}")
        else:
            print("- None")
        print("") # Add a newline for better readability

# Execute the classification function
classify_fusion_experiments()

# The final answer in the requested format will summarize the classification.
# For each category, I will list the numbers corresponding to the experiments.
# 1) Tokamaks -> A, B, and C
# 2) LHD -> B and C
# 3) Wendelstein 7-X -> B and C
# 4) NCSX -> B and C
# 5) Reversed Field Pinches -> A, B, and C
# This can't be represented as a single character or number.
# Let's format it as a string.
final_answer = """
Only B: None
B and C: LHD, Wendelstein 7-X, NCSX
A, B, and C: Tokamaks, Reversed Field Pinches
"""
# The format seems to expect a single value, I will output the most complete group's name.
# Or I could list the classification for each number. Let's try that.
final_answer_formatted = "1: (A,B,C), 2: (B,C), 3: (B,C), 4: (B,C), 5: (A,B,C)"
# Since I must choose one response, I will provide the full text block answer.

# However, the prompt asks me to "directly return the answer with the format <<<answer content>>>".
# Given the complexity, a string seems most appropriate.

print("<<<Only B: None; B and C: LHD, Wendelstein 7-X, NCSX; A, B, and C: Tokamaks, Reversed Field Pinches>>>")
