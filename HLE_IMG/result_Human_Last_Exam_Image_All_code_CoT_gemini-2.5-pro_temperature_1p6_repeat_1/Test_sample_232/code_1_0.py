def solve_medical_puzzle():
    """
    This function identifies the atrophied organ part and the corresponding disorder
    based on the provided medical image.
    """
    # The organ part that has undergone significant atrophy
    atrophied_part = "caudate nucleus"
    
    # The name of the disorder associated with this type of atrophy
    disorder_name = "Huntington's disease"
    
    # Combine the two parts into the required format "organ, disorder"
    final_answer = f"{atrophied_part}, {disorder_name}"
    
    print(final_answer)

solve_medical_puzzle()