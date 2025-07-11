def solve_medical_puzzle():
    """
    This function identifies the atrophied organ part and the associated disorder from the image.
    The image shows marked atrophy of the caudate nucleus, a classic sign of Huntington's disease.
    """
    
    # The part of the organ that has atrophied
    atrophied_part = "caudate nucleus"
    
    # The name of the disorder
    disorder = "Huntington's disease"
    
    # Format the answer as "organ, disorder"
    # Note: "caudate nucleus" is the specific atrophied part of the brain.
    final_answer = f"{atrophied_part}, {disorder}"
    
    print(final_answer)

solve_medical_puzzle()