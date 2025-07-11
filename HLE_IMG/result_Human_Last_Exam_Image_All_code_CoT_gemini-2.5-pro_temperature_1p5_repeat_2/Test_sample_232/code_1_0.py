def identify_disorder_from_image():
    """
    Identifies the atrophied organ part and the associated disorder based on the provided image.
    The image shows classic signs of hippocampal atrophy.
    """
    
    # The part of the organ that has atrophied
    atrophied_part = "hippocampus"
    
    # The disorder associated with this atrophy
    disorder = "Alzheimer's disease"
    
    # Print the result in the format "organ, disorder"
    print(f"{atrophied_part}, {disorder}")

identify_disorder_from_image()