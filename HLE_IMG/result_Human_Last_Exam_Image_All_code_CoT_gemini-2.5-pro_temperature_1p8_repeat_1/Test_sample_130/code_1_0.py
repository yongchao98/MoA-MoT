import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
string_io = io.StringIO()
# Redirect stdout to the string buffer
sys.stdout = string_io

def identify_element_from_spectrum():
    """
    Identifies an element based on the visual characteristics of its emission spectrum.
    """
    element_name = "Iron (Fe)"
    reasoning = [
        "1. Each element possesses a unique emission spectrum, which acts as its 'fingerprint'.",
        "2. The provided spectrum is an emission spectrum, characterized by bright lines on a dark background.",
        "3. The spectrum in the image is extremely complex, containing a very large number of emission lines across the entire visible range (from blue to red).",
        "4. Simple elements like Hydrogen or Helium have much simpler spectra with fewer lines.",
        "5. This high density of lines is a hallmark of heavier elements, particularly transition metals.",
        f"6. By comparing this spectral pattern to known reference spectra, it is identified as that of {element_name}."
    ]

    print(f"The element with the spectral lines shown in the image is: {element_name}\n")
    print("Reasoning:")
    for point in reasoning:
        print(point)

identify_element_from_spectrum()

# Get the content from the string buffer
output = string_io.getvalue()
# Restore original stdout
sys.stdout = original_stdout

# Print the captured output
print(output)
print("<<<Iron>>>")