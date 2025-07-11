import sys

def solve():
    """
    This script identifies the atrophied brain structure and the corresponding disorder from the provided image.
    The image shows a coronal section of a brain with marked atrophy of the caudate nucleus,
    leading to the characteristic 'boxcar' appearance of the lateral ventricles.
    This is a hallmark sign of Huntington's disease.
    """
    
    # The part of the organ that has atrophied
    atrophied_structure = "caudate nucleus"
    
    # The name of the disorder
    disorder = "Huntington's disease"
    
    # Combine them into the required format "organ, disorder"
    final_answer = f"{atrophied_structure}, {disorder}"
    
    print(final_answer)

solve()