def get_piano_fingering():
  """
  Provides the most comfortable fingering for the given piano chord.
  """
  # The chord notes and their corresponding finger numbers (Right Hand)
  # 1=Thumb, 2=Index, 3=Middle, 4=Ring, 5=Pinky
  notes = ["Bb3", "Db4", "F4", "Bb4", "Db5"]
  fingers = [1, 2, 3, 4, 5]

  # Format the output string as "Note Finger, Note Finger, ..."
  # The f-string builds the final string with each note and its finger.
  fingering_chart = f"{notes[0]} {fingers[0]}, {notes[1]} {fingers[1]}, {notes[2]} {fingers[2]}, {notes[3]} {fingers[3]}, {notes[4]} {fingers[4]}"
  
  print(fingering_chart)

get_piano_fingering()