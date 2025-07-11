import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a new StringIO object
sys.stdout = captured_output = io.StringIO()

def identify_structure():
    """
    Identifies a structure from satellite imagery based on geographic coordinates.
    """
    # 1. Define the coordinates in Degrees, Minutes, Seconds (DMS)
    lat_d, lat_m, lat_s = 29, 6, 18.75
    lon_d, lon_m, lon_s = 103, 47, 50.28

    # 2. Convert DMS to Decimal Degrees (DD)
    # The formula is: DD = Degrees + (Minutes/60) + (Seconds/3600)
    # For South and West coordinates, the decimal value is negative.
    latitude_dd = lat_d + (lat_m / 60) + (lat_s / 3600)
    longitude_dd = -1 * (lon_d + (lon_m / 60) + (lon_s / 3600))

    # 3. Information about the site based on research at these coordinates
    site_name = "Presidio de San Vicente"
    other_names = "El Polvo or the Scholes Site"
    description = "The ruins of an 18th-century Spanish fort."
    history = "Built around 1773-1774, it was part of a line of military forts (presidios) established by the Spanish Empire along the northern frontier of New Spain to defend against raids."
    is_historic_landmark = "Yes, it is a significant historic and archaeological site in Presidio County, Texas."

    # 4. Print the findings
    print("--- Analysis of the Structure at 29°06'18.75''N, 103°47'50.28''W ---")
    print("\n[Step 1 & 2: Coordinate Conversion]")
    print("The provided coordinates are converted from DMS to Decimal Degrees for analysis.")
    print("Latitude conversion equation:")
    print(f"{lat_d} + ({lat_m} / 60) + ({lat_s} / 3600) = {latitude_dd:.6f}°")
    print("Longitude conversion equation:")
    print(f"-1 * ({lon_d} + ({lon_m} / 60) + ({lon_s} / 3600)) = {longitude_dd:.6f}°")

    print("\n[Step 3 & 4: Site Identification]")
    print("\nWhat is the square ancient structure?")
    print(f"The structure is the {description}")
    print(f"It is known as '{site_name}', and also referred to as '{other_names}'.")

    print("\nIs it a historic landmark?")
    print(f"{is_historic_landmark}")
    print(f"Historical Context: {history}")

identify_structure()

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
print(output)

# Final answer in the required format
final_answer = "The structure is the ruins of the 18th-century Spanish fort, Presidio de San Vicente. It is a historic landmark."
print(f"<<<{final_answer}>>>")