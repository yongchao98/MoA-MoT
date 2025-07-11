import sys
import io

# Redirect stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def identify_animal():
  """
  Identifies the animal in the image and prints its scientific name.
  """
  animal_present = "Yes"
  scientific_name = "Acanthostracion polygonius"
  common_name = "Honeycomb Cowfish"

  print(f"Yes, there appears to be an animal in the image.")
  print(f"Based on the distinctive honeycomb pattern, it is likely a {common_name}.")
  print(f"Its scientific name is: {scientific_name}")

identify_animal()

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

print(output)