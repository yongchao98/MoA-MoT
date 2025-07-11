def solve_fingering():
  """
  Determines and prints the most comfortable right-hand piano fingering for the chord:
  Bb3, Db4, F4, Bb4, Db5
  
  The standard piano fingering numbers for the right hand are:
  1 = Thumb
  2 = Index Finger
  3 = Middle Finger
  4 = Ring Finger
  5 = Pinky Finger
  """
  
  # The chord requires a wide stretch from Bb3 to Db5.
  # The most comfortable and standard fingering uses all five fingers in order.
  
  # Note: Bb3, Finger: 1 (Thumb)
  note1 = "Bb3"
  finger1 = 1
  
  # Note: Db4, Finger: 2 (Index)
  note2 = "Db4"
  finger2 = 2
  
  # Note: F4, Finger: 3 (Middle)
  note3 = "F4"
  finger3 = 3

  # Note: Bb4, Finger: 4 (Ring)
  note4 = "Bb4"
  finger4 = 4
  
  # Note: Db5, Finger: 5 (Pinky)
  note5 = "Db5"
  finger5 = 5
  
  # The final answer is formatted as a single string.
  final_answer = f"{note1} {finger1}, {note2} {finger2}, {note3} {finger3}, {note4} {finger4}, {note5} {finger5}"
  print(final_answer)

solve_fingering()