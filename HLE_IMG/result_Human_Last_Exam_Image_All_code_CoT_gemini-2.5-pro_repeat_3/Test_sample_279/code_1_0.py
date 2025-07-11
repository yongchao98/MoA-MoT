def get_caries_classification(location_description):
    """
    Determines the G.V. Black caries classification based on a description.
    """
    classification_map = {
        "pits_and_fissures_on_occlusal_surface": "Class I",
        "proximal_surface_of_posterior_tooth": "Class II",
        "proximal_surface_of_anterior_tooth": "Class III",
        "proximal_surface_of_anterior_tooth_involving_incisal_edge": "Class IV",
        "cervical_third_of_tooth": "Class V",
        "cusp_tips_or_incisal_edges": "Class VI"
    }
    
    return classification_map.get(location_description, "Unknown")

# Step 1: Analyze the image to determine the lesion's location.
# The image shows decay (and a filling) on the chewing surface of a molar.
# This area is known as the occlusal surface, and the decay is in the pits and fissures.
lesion_location = "pits_and_fissures_on_occlusal_surface"

# Step 2: Get the classification code.
classification_code = get_caries_classification(lesion_location)

# Step 3: Print the result.
print("According to the G.V. Black classification system:")
print(f"Lesion Location: On the pits and fissures of the occlusal surface.")
print(f"This corresponds to a Classification Code: {classification_code}")
print("\nExplanation: Class I cavities are defined as those occurring in the pits and fissures of the teeth, primarily on the occlusal surfaces of molars and premolars. The image shows an existing restoration with recurrent decay in this exact location.")
