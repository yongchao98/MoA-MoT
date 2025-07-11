import sys
import io

# Redirect stdout to a string buffer
old_stdout = sys.stdout
sys.stdout = mystdout = io.StringIO()

def identify_animal():
  """
  Identifies the animal in the image and prints its scientific name.
  """
  animal_present = "Yes, there is an animal in this image."
  animal_type = "It is an insect, specifically a type of longhorn beetle."
  scientific_name = "Stictoleptura cordigera"

  print(animal_present)
  print(f"{animal_type}")
  print(f"The scientific name of this beetle is: {scientific_name}")

identify_animal()

# Restore stdout and get the output
sys.stdout = old_stdout
output = mystdout.getvalue()

# The user has requested not to be asked to copy and paste.
# The following code simulates the final output they would see.
# In a real execution, only the print statements within the function would run.
# I am printing the captured output to fulfill the instructions.
print(output)