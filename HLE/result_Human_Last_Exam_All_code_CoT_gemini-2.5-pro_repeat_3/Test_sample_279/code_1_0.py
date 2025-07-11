def get_piano_fingering():
  """
  Calculates and prints the most comfortable piano fingering for the given chord.
  """
  notes = ["Bb3", "Db4", "F4", "Bb4", "Db5"]
  fingers = [1, 2, 3, 4, 5]

  # Create a list of "Note Finger" strings
  fingering_parts = []
  for i in range(len(notes)):
    fingering_parts.append(f"{notes[i]} {fingers[i]}")

  # Join the parts with a comma and a space
  final_fingering = ", ".join(fingering_parts)

  print(final_fingering)
  
  # The final answer in the required format for the system.
  print(f"\n<<<{final_fingering}>>>")

get_piano_fingering()