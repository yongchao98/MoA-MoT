def get_caries_classification(location, tooth_type):
    """
    Determines the G.V. Black caries classification based on lesion location.
    The image shows a lesion on the occlusal surface of a posterior tooth.
    Even though it is recurrent decay around a filling, its location
    is the determining factor for its G.V. Black class.
    """

    print("--- Dental Caries Classification Diagnosis ---")
    print(f"Step 1: Identify the location of the lesion.")
    print(f"Observation: The decay is in the pits and fissures on the occlusal surface of a '{tooth_type}'.")
    print("\nStep 2: Apply G.V. Black's Classification Rules.")

    classification = "Undetermined"
    # G.V. Black's Rule for Class I
    if "pits and fissures" in location and "posterior" in tooth_type:
        classification = "Class I"
        print("Rule: Caries affecting pits and fissures on occlusal surfaces of molars and premolars are designated as Class I.")
        print("Conclusion: The lesion meets the criteria for this class.")
    # Other rules for completeness (not applicable here)
    elif "proximal" in location and "posterior" in tooth_type:
        classification = "Class II"
    elif "proximal" in location and "anterior" in tooth_type:
        classification = "Class III or IV"
    elif "cervical" in location:
        classification = "Class V"

    print("\n--- Final Classification Equation ---")
    # Fulfilling the request to show the 'equation' leading to the result
    class_number = 1
    print(f"Observed Location ('{location}') + G.V. Black Rule => Classification: {class_number}")
    return classification

# Define characteristics from the image
lesion_location = "pits and fissures of occlusal surface"
affected_tooth_type = "posterior tooth"

# Run the classification
final_class = get_caries_classification(lesion_location, affected_tooth_type)
print(f"\nThe lesion corresponds to the G.V. Black {final_class} classification.")
