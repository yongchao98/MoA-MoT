def find_newtonian_note():
    """
    Solves for the musical note based on Newton's color circle and a derived equation.
    """

    # Step 1: Define the mappings from color to note and note to a numerical scale position (C-Major).
    color_note_map = {
        'Yellow': 'E',
        'Blue': 'G',
        'Green': 'F' # The color produced by mixing Yellow and Blue
    }

    note_number_map = {
        'E': 3,
        'F': 4,
        'G': 5
    }

    # Step 2: Identify the colors and corresponding notes from the problem.
    # The start of the range is the color made by mixing Yellow and Blue.
    start_color_component1 = 'Yellow'
    start_color_component2 = 'Blue'
    start_color_mixed = 'Green'

    # Retrieve the note names and their numerical values.
    note1_name = color_note_map[start_color_component1]
    note2_name = color_note_map[start_color_component2]
    result_note_name = color_note_map[start_color_mixed]

    note1_num = note_number_map[note1_name]
    note2_num = note_number_map[note2_name]
    result_note_num = note_number_map[result_note_name]

    # Step 3: Explain the logic and print the equation with all its numbers.
    # The note for the mixed color 'Green' corresponds to the start of the specified range.
    print(f"The specified range starts with the color made from '{start_color_component1}' and '{start_color_component2}', which is '{start_color_mixed}'.")
    print(f"On Newton's circle, these colors correspond to notes E ({note1_num}), G ({note2_num}), and F ({result_note_num}).")
    print("\nThe note for Green (F) can be found by averaging the numerical positions of the notes for Yellow (E) and Blue (G).")
    print("\nFinal Equation:")
    
    # Print the equation showing each number as requested.
    # The first number is 3 (for Yellow/E).
    # The second number is 5 (for Blue/G).
    # The result of the division is 4 (for Green/F).
    print(f"({note1_num} + {note2_num}) / 2 = {int((note1_num + note2_num) / 2)}")
    
    print(f"\nThe note corresponding to the number {result_note_num} is {result_note_name}.")

# Execute the function to find and print the answer.
find_newtonian_note()