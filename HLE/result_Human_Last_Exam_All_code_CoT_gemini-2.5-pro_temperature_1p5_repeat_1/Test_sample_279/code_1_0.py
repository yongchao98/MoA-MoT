def get_piano_fingering():
  """
  Determines and prints the most comfortable right-hand fingering for the specified piano chord.
  
  The chord notes are: Bb3, Db4, F4, Bb4, Db5.
  The fingering is determined based on standard piano technique principles for wide chords.
  """
  
  # Notes of the chord
  notes = ["Bb3", "Db4", "F4", "Bb4", "Db5"]
  
  # Standard fingering for a wide chord like this (1=thumb, 5=pinky)
  fingers = [1, 2, 3, 4, 5]
  
  # Combine notes and fingerings into the desired format
  # e.g., "Bb3 1, Db4 2, F4 3, Bb4 4, Db5 5"
  output_parts = []
  for i in range(len(notes)):
    part = f"{notes[i]} {fingers[i]}"
    output_parts.append(part)
    
  final_string = ", ".join(output_parts)
  
  print(final_string)

get_piano_fingering()