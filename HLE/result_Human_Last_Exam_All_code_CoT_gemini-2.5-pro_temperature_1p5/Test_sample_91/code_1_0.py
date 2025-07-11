# The scale identified in Erroll Garner's right-hand melody
# between 0:39 and 0:43 of "All My Loves Are You" is the C Natural Minor scale.

# This scale is composed of the following notes.
# The code below defines and prints these notes.

def get_scale_notes(root_note, scale_type):
    """
    This function returns the notes of a specific scale.
    For this problem, we will focus on C Natural Minor.
    """
    if root_note == "C" and scale_type == "Natural Minor":
        # The notes of the C Natural Minor scale are C, D, Eb, F, G, Ab, Bb
        return ["C", "D", "Eb", "F", "G", "Ab", "Bb"]
    else:
        return []

def main():
    """
    Main function to identify and print the scale notes.
    """
    scale_name = "C Natural Minor"
    root = "C"
    scale_type = "Natural Minor"

    print(f"The scale played by Erroll Garner in the specified passage is the {scale_name} scale.")
    print("The notes that make up this scale are:")

    notes = get_scale_notes(root, scale_type)
    
    # Print each note in the scale, fulfilling the prompt's requirement
    # to output each component part.
    print(" + ".join(notes))

if __name__ == "__main__":
    main()
