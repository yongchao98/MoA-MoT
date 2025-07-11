def get_piano_fingering():
  """
  Calculates and prints the most comfortable right-hand fingering for a given chord.
  """
  notes = ["Bb3", "Db4", "F4", "Bb4", "Db5"]
  # Based on ergonomic principles for a wide-spanned chord with black and white keys,
  # the most direct fingering is the most comfortable, assuming the hand can make the stretch.
  fingers = ["1", "2", "3", "4", "5"]

  # Assemble the formatted string
  output_parts = []
  for i in range(len(notes)):
    note = notes[i]
    finger = fingers[i]
    # Each part of the output includes the note and its corresponding finger
    output_parts.append(f"{note} {finger}")
  
  final_answer = ", ".join(output_parts)
  print(final_answer)

get_piano_fingering()