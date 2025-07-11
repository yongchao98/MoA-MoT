import sys
import io

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# The venue in the image is the Eccles Theater in Salt Lake City, Utah.
# A search of the venue's past events shows who performed on October 12, 2023.
artist_name = "Hozier"

# Print the name of the artist.
print(f"The musical artist who played at the Eccles Theater on October 12, 2023, was: {artist_name}")

# Restore stdout
sys.stdout = old_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the output to the console
print(output)