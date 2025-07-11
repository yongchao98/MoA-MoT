import sys

def solve():
    """
    This function identifies the atrophied structure and the disorder shown in the medical image.
    """
    # The image shows a coronal section of a brain.
    # There is marked atrophy of the caudate nucleus, which is a part of the basal ganglia.
    # This specific atrophy leads to the enlargement of the lateral ventricles, a hallmark sign.
    atrophied_structure = "Caudate nucleus"
    
    # This pathological finding is characteristic of a specific neurodegenerative disorder.
    disorder = "Huntington's disease"
    
    # Printing the answer in the format "atrophied part, disorder".
    print(f"{atrophied_structure}, {disorder}")

solve()