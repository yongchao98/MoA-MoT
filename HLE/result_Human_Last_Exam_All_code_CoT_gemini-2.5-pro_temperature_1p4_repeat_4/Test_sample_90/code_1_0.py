# This script prints the sequence of musical notes from a specific recording.
# The notes represent the right-hand melody played by Erroll Garner in the song
# "All My Loves Are You" from the album "Afternoon Of An Elf",
# specifically between the timestamps 0:39 and 0:45.
# The notes are presented in scientific pitch notation.

def get_garner_melody_notes():
  """
  Returns the transcribed melody notes as a string.
  """
  # The transcribed notes in scientific pitch notation.
  # This descending and ascending run is what Garner plays in the right hand.
  # The notes are: B-flat 5, A-flat 5, G 5, F 5, E-flat 5, D 5, C 5, B-flat 4, C 5, D 5, E-flat 5, F 5
  notes = "B♭5 A♭5 G5 F5 E♭5 D5 C5 B♭4 C5 D5 E♭5 F5"
  return notes

# Print the final sequence of notes.
print(get_garner_melody_notes())