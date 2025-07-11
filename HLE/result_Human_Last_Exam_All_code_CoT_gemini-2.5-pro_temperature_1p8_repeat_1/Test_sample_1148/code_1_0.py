def find_borges_names():
  """
  Identifies the three individuals mentioned by Jorge Luis Borges in his essay
  "Our Poor Individualism" to illustrate national patriotism.
  """
  # The Greek example refers to the dialogue "Crito," written by Plato,
  # where Socrates argues he must obey the laws of Athens.
  greek = "Plato"

  # The English example is the famous lexicographer and writer
  # credited with the quote: "Patriotism is the last refuge of a scoundrel."
  englishman = "Samuel Johnson"

  # The German example is the philosopher who gave the "Addresses to the German Nation."
  german = "Johannes Fichte"

  # Combine the names into a single string, separated by commas.
  result = f"{greek}, {englishman}, {german}"

  print(result)

find_borges_names()