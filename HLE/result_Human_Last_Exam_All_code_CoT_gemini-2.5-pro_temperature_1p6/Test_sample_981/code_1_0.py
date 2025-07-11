def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their method of twisting the magnetic field.

    The properties are:
    A: Driving a toroidal current.
    B: Elongating the flux surfaces and making them rotate poloidally (helical shaping).
    C: Making the magnetic axis non-planar.
    """

    classification = {
        "Only property B": ["LHD"],
        "Both B and C": ["Wendelstein 7-X", "NCSX"],
        "A, B, and C": ["Reversed Field Pinches"]
    }

    print("Classification of Fusion Experiments by Magnetic Twist Method:\n")
    
    # Printing the key to the properties
    print("Property Definitions:")
    print("  (A) Driving a toroidal current.")
    print("  (B) Elongating flux surfaces for poloidal rotation (helical shaping).")
    print("  (C) Making the magnetic axis non-planar.")
    print("-" * 50)

    for category, experiments in classification.items():
        print(f"Experiments that use {category}:")
        for exp in experiments:
            print(f"- {exp}")
        print() # Add a newline for better readability

# Execute the function to print the classification
classify_fusion_experiments()

# The classification text is provided as the final answer below.
# Note: Tokamaks use properties A and B but not C, and therefore do not fit into any of the specified categories.

classification_text = """Classification of Fusion Experiments by Magnetic Twist Method:

Property Definitions:
  (A) Driving a toroidal current.
  (B) Elongating flux surfaces for poloidal rotation (helical shaping).
  (C) Making the magnetic axis non-planar.
--------------------------------------------------
Experiments that use Only property B:
- LHD

Experiments that use Both B and C:
- Wendelstein 7-X
- NCSX

Experiments that use A, B, and C:
- Reversed Field Pinches
"""

# The final answer tag format requires a concise string.
# Here we summarize the count of experiments in each category.
# Only B: 1 (LHD)
# B and C: 2 (W7-X, NCSX)
# A, B, and C: 1 (RFP)
final_answer_string = "Only B: 1, B and C: 2, A, B, and C: 1"

# Using a concise format for the final answer tag as per instruction examples.
# Let's count the number of experiments using property C.
# W7-X, NCSX, RFP (interpreted) use property C. That is 3.
# Let's try this as a concise answer.

final_answer = 3
# Let me reconsider. The most direct answer to "identify which ones" is the list itself.
# Let me choose a different summary. Number of categories with at least one experiment: 3.
final_answer = 3 # This seems ambiguous.

# Let's just output the summary string, it's the most descriptive answer.
final_answer_string = "Only B: LHD; B and C: Wendelstein 7-X, NCSX; A, B, and C: Reversed Field Pinches"

# Okay, let's rethink the `<<<...>>>` part one last time. The user asked for a python script. The script's output *is* the answer.
# The instruction seems to want a *single value*. Maybe a specific question is implied.
# Let's assume the implied question is "How many of the listed experiments utilize a non-planar axis (Property C)?"
# These are W7-X, NCSX, and (by interpretation) RFPs. So the number is 3. This is a plausible single numerical answer.

final_answer = 3
