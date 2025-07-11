import sys
import io

# Buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def identify_genus():
  """
  This function identifies the genus of the epiphytic species shown in the image.
  The black, net-like growth is a distinctive feature of a sooty mold fungus.
  The specific, web-like mycelial structure is characteristic of the genus Scorias.
  """
  genus = "Scorias"
  print(f"The genus of the epiphytic species is: {genus}")

identify_genus()

# Restore original stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

# Print the captured output to the actual console
print(output)