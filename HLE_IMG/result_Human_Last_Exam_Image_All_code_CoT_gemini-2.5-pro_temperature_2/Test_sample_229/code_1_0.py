def get_song_info():
  """
  This function provides information about the song shown in the image.
  
  The analysis is based on the following observations from the sheet music:
  - Time Signature: 4/4
  - Key/Chords: Based in C, uses standard jazz chords (C, E7, F, G7, Dm, etc.), suggesting a jazz or popular music style.
  - Melody & Rhythm: The melody is highly syncopated, which is a key characteristic of Ragtime and Dixieland jazz.
  - Title: Faint text at the top appears to read "Muskrat...", pointing to the well-known standard "Muskrat Ramble".
  - Confirmation: Comparing the melody and chords with known lead sheets for "Muskrat Ramble" confirms the identification. This tune was composed by Kid Ory in 1926.
  - Origin & Style: Kid Ory was a foundational figure in New Orleans jazz, placing the song's origin in the USA and its style as Dixieland.
  """
  title = "Muskrat Ramble"
  musical_style = "Dixieland Jazz"
  country_of_origin = "United States"

  print(f"Title: {title}")
  print(f"Musical Style: {musical_style}")
  print(f"Country of Origin: {country_of_origin}")

get_song_info()