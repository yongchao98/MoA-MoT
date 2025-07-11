def identify_dmPFC_subregion():
    """
    Identifies the dmPFC sub-region most purely activated by emotion processing
    by analyzing the primary functions of each parcellation based on neuroscience literature.
    """
    
    # Step 1: Represent the functional specialization of the four dmPFC sub-regions.
    # This data is based on large-scale meta-analyses of fMRI studies.
    # The first function listed is considered the most primary or defining role.
    functional_parcellations = {
        "caudal-right": {
            "choice": "A",
            "primary_function": "Motor/Cognitive Control",
            "description": "Primarily involved in action selection, response inhibition, and monitoring outcomes. Its activation during emotional tasks is often related to preparing or executing a response."
        },
        "rostroventral": {
            "choice": "B",
            "primary_function": "Emotion and Valuation",
            "description": "This region is core to generating and regulating emotional states, assessing the value of stimuli, and affective self-reference. Its function is most selectively tied to emotion."
        },
        "rostrodorsal": {
            "choice": "C",
            "primary_function": "Social Cognition (Mentalizing)",
            "description": "Specialized for thinking about the mental states of others (Theory of Mind) and complex, cognitive self-reflection. While related to social emotions, its core role is more cognitive."
        },
        "caudal-left": {
            "choice": "D",
            "primary_function": "Motor/Cognitive Control",
            "description": "Similar to its right-hemisphere counterpart, this area is primarily involved in action and cognitive control, not pure emotional processing."
        }
    }
    
    print("Evaluating dmPFC sub-regions for purity of emotion processing...\n")
    
    # Step 2 & 3: Define the criterion and analyze each region.
    # The criterion is "most purely activate by emotion", meaning emotion is the primary function.
    purest_region_name = None
    purest_region_choice = None
    
    for region, data in functional_parcellations.items():
        print(f"Region: {region} (Choice: {data['choice']})")
        print(f"  - Primary Function: {data['primary_function']}")
        print(f"  - Analysis: {data['description']}")
        
        # Check if the primary function is 'Emotion and Valuation'. This signifies purity.
        if data["primary_function"] == "Emotion and Valuation":
            purest_region_name = region
            purest_region_choice = data['choice']
        print("-" * 60)
        
    # Step 4: Output the conclusion.
    print("\nConclusion:")
    print("The caudal regions are primarily for motor and cognitive control.")
    print("The rostrodorsal region is primarily for cognitive social processing (mentalizing).")
    print("The rostroventral region is the only sub-region whose primary, defining role is directly related to emotion and value processing.")
    print(f"\nTherefore, the sub-region most purely activated by emotion is '{purest_region_name}'.")
    print(f"The corresponding answer choice is: {purest_region_choice}")

# Run the analysis
identify_dmPFC_subregion()