import sys
import io

# Redirect stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def identify_painting():
  """
  Identifies and prints the details of the provided painting.
  The image shows a dramatic scene with a windmill under a stormy sky,
  and a fire in the background. This scene is characteristic of a specific
  historical painting from the 19th century.
  """
  artist = "Théodore Gudin"
  title = ("The Burning of the 'Faubourg de l'Esplan' Windmill, "
           "during the siege of the Citadel of Antwerp, 4 December 1832")
  year_of_event = 1832

  print(f"Artist: {artist}")
  print(f"Title: {title}")
  print(f"The event depicted occurred in the year: {year_of_event}")

identify_painting()

# Restore stdout and get the captured content
sys.stdout = old_stdout
output = captured_output.getvalue()

# The final response should still present the information clearly.
# Let's print the captured output for the user.
print(output)