def find_musical_joke_measures():
  """
  Identifies and prints the details of a famous musical joke
  in a piano concerto by Saint-Saëns.
  """
  
  # Define the components of the answer
  composer_surname = "Saint-Saëns"
  opus_number = 22
  start_measure = 113
  end_measure = 116
  
  # Format the output string as requested: "composer, opus, start-end"
  # The flourish begins in measure 113 and the recovery completes in measure 116.
  output = f"{composer_surname}, {opus_number}, {start_measure}-{end_measure}"
  
  print(output)

find_musical_joke_measures()