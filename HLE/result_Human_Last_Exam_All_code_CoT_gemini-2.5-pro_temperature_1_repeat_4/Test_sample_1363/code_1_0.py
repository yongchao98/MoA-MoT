def find_the_dance():
    """
    This function analyzes the technical principles of standard ballroom dances
    to answer a specific question about turning technique.
    """
    question = "In which dance is it impossible to overturn a reverse turn without disregarding the technique?"
    
    choices = {
        'A': 'Viennese Waltz',
        'B': 'English Waltz',
        'C': 'European Tango',
        'D': 'Slow Foxtrot',
        'E': 'Quickstep'
    }

    # Analysis of the dance techniques:
    # The dances listed can be categorized by their core movement principles.
    # Waltz, Viennese Waltz, Foxtrot, and Quickstep are "swing" dances. They utilize
    # techniques like Contra Body Movement (CBM), sway, and rise & fall to generate
    # flowing, continuous, and rotational movement. While overturning a turn in these
    # dances is technically incorrect, it is achieved by exaggerating or misusing
    # these inherent swing mechanics.
    #
    # European Tango is fundamentally different. Its technique is characterized by
    # staccato action, no rise & fall, and no sway. Movement is sharp, direct, and
    # kept under the body's center. Turns are executed with precise footwork and
    # sharp head turns, not with a continuous, swinging rotation of the body.
    #
    # To "overturn" a turn in Tango, one would need to introduce a rotational swing or sway,
    # which are not part of the Tango's technical foundation. Doing so is not merely
    # poor technique; it is a complete disregard for the defining characteristics of the dance.
    
    correct_answer_key = 'C'
    explanation = (
        "In European Tango, the technique forbids sway and the continuous rotational "
        "body flight found in 'swing' dances. A reverse turn is executed with a sharp, "
        "contained action. Attempting to 'overturn' it would necessitate introducing "
        "a swinging motion and sway, which are fundamentally contrary to the core "
        "principles of Tango. Therefore, it's considered impossible to do so without "
        "completely disregarding the dance's technique."
    )

    print("The question is:")
    print(f"'{question}'")
    print("\nBased on an analysis of ballroom dance techniques:")
    print(explanation)
    print(f"\nThe correct option is: {correct_answer_key} - {choices[correct_answer_key]}")

find_the_dance()