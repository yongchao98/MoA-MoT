def find_languages():
  """
  This function identifies two languages based on a set of linguistic clues
  and prints the names of the languages.
  """

  # Language a: No 'k' or 'w', has 'à', is a recognized minority language.
  # This points to Scottish Gaelic.
  language_a = "Scottish Gaelic"

  # Language b: Uses "ggj" and "skt", is a recognized minority language.
  # The "ggj" cluster is characteristic of Northern Sami.
  language_b = "Northern Sami"

  print(f"Language a is: {language_a}")
  print(f"Language b is: {language_b}")

find_languages()