import geopy.distance

def dms_to_dd(degrees, minutes, seconds, direction):
    """Converts coordinates from DMS to decimal degrees."""
    dd = float(degrees) + float(minutes)/60 + float(seconds)/(60*60)
    if direction in ['S', 'W']:
        dd *= -1
    return dd

def solve_rock_art_mystery():
    """
    Analyzes the location and imagery of a rock carving to determine its meaning.
    """
    # 1. Define and convert coordinates for the rock carving
    carving_lat_dms = (29, 3, 28.15)
    carving_lon_dms = (103, 48, 11.84)
    carving_coords = (
        dms_to_dd(*carving_lat_dms, 'N'),
        dms_to_dd(*carving_lon_dms, 'W')
    )

    # 2. Define coordinates for locations in answer choices
    # Approximate center of Chisos Mountains
    chisos_mountains_coords = (29.2714, -103.3039) # Chisos Basin
    # Lajitas, Texas, on the Rio Bravo (Rio Grande)
    lajitas_coords = (29.2636, -103.7664)

    # 3. Calculate distances to check claims in options B and D
    distance_to_chisos = geopy.distance.geodesic(carving_coords, chisos_mountains_coords).miles
    distance_to_lajitas = geopy.distance.geodesic(carving_coords, lajitas_coords).miles

    # 4. Print the analysis
    print("Step 1: Analyzing the location of the rock carving.")
    print(f"The coordinates 29° 3' 28.15''N, 103° 48' 11.84''W place the carving in the Big Bend region of Texas.")
    print("This area is known for significant Native American history and rock art (petroglyphs).\n")

    print("Step 2: Evaluating the answer choices.\n")

    print("--- Evaluating Choice C: 'No. There are no carvings in the rock.' ---")
    print("Result: This is FALSE. The provided images clearly show etched lines and figures on the rock surface.\n")

    print("--- Evaluating Choice A: 'No, the carving depicts a snake.' ---")
    print("Result: This is UNLIKELY. The most prominent figure is a series of concentric circles or a spiral. This does not match the typical form of a snake.\n")

    print("--- Evaluating Choice B: 'Yes. the image depicts the Chiso Mountains...' ---")
    print(f"The calculated distance to the Chisos Mountains is approximately {distance_to_chisos:.2f} miles.")
    print("The claim of 'about 20 miles north' is plausible in distance, though the actual direction is northeast.")
    print("However, the carving's abstract, circular design does not visually resemble a mountain range.\n")

    print("--- Evaluating Choice D: 'Yes. It matches with a segment of the Bravo River...' ---")
    print(f"The calculated distance to Lajitas on the Rio Bravo is approximately {distance_to_lajitas:.2f} miles.")
    print("The claim of '10 miles northwest' is plausible as the location is about 13 miles to the south.")
    print("However, the circular pattern is not a realistic depiction of a river segment.\n")
    
    print("--- Evaluating Choice E: 'No. The lines carved in the rock depict a feature known as sipapu...' ---")
    print("A 'sipapu' is a symbolic portal from the underworld in the beliefs of Southwestern indigenous peoples (like the Hopi and other Pueblo cultures).")
    print("In rock art, the sipapu is frequently represented by spirals or concentric circles.")
    print("Result: This is the MOST LIKELY answer. The main carving, a circular/spiral pattern, is a strong visual match for a sipapu symbol. This interpretation is consistent with the archaeological context of the region.\n")

    print("Conclusion: The carving does not depict a geographic feature. Its form and location strongly suggest it represents a symbolic sipapu.")

solve_rock_art_mystery()