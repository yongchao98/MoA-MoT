def identify_garner_scale():
    """
    Identifies and describes the scale used by Erroll Garner in a specific
    passage of "All My Loves Are You".
    """
    scale_name = "C Lydian"
    
    # The notes of the C major scale for reference
    c_major_scale = ['C', 'D', 'E', 'F', 'G', 'A', 'B']
    
    # The notes Erroll Garner uses in the run between 0:39-0:43
    # When arranged in ascending order, they form the C Lydian scale.
    c_lydian_scale = ['C', 'D', 'E', 'F#', 'G', 'A', 'B']

    print("Analysis of the piano melody from 0:39 to 0:43:")
    print(f"The passage is played over a C major chord.")
    print(f"The notes used in the fast descending run are characteristic of the {scale_name} scale.")
    print("This scale is like a major scale but with a raised fourth note, giving it a bright, distinctive sound common in jazz.")
    
    # Fulfilling the requirement to output each 'number' (note) in the 'equation' (scale).
    print("\nThe scale's notes are:")
    final_output = " -> ".join(c_lydian_scale)
    print(final_output)

identify_garner_scale()