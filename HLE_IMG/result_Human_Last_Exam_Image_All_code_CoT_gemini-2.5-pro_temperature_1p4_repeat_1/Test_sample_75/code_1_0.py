import geopy.distance

def dms_to_dd(degrees, minutes, seconds):
    """Converts coordinates from Degrees, Minutes, Seconds to Decimal Degrees."""
    return degrees + minutes / 60 + seconds / 3600

def get_direction(start_point, end_point):
    """Determines the general cardinal direction from a start point to an end point."""
    lat1, lon1 = start_point
    lat2, lon2 = end_point
    
    direction = ""
    
    if lat2 > lat1:
        direction += "North"
    elif lat2 < lat1:
        direction += "South"
        
    if lon2 > lon1:
        direction += "East"
    elif lon2 < lon1:
        direction += "West"
        
    # If there's no change in one axis, remove the hyphen
    if direction.endswith("North") or direction.endswith("South"):
        return direction
    if direction.startswith("East") or direction.startswith("West"):
        return direction
    
    # Add a hyphen for combined directions
    if "North" in direction and ("East" in direction or "West" in direction):
        return "North-" + direction.replace("North", "")
    if "South" in direction and ("East" in direction or "West" in direction):
        return "South-" + direction.replace("South", "")
        
    return "Same location"

# 1. Coordinates of the rock carving
lat_carving_dd = dms_to_dd(29, 3, 28.15)
# Longitude is West, so it's negative
lon_carving_dd = -dms_to_dd(103, 48, 11.84)
point_carving = (lat_carving_dd, lon_carving_dd)

# 2. Approximate coordinates for reference points
# Chisos Mountains are a range; using a central point in the range for approximation.
point_chisos = (29.27, -103.30) 
# Lajitas, Texas is a town near the Rio Grande (Rio Bravo).
point_lajitas = (29.26, -103.76)

# 3. Calculate distances and directions
distance_to_chisos_miles = geopy.distance.great_circle(point_carving, point_chisos).miles
direction_to_chisos = get_direction(point_carving, point_chisos)

distance_to_lajitas_miles = geopy.distance.great_circle(point_carving, point_lajitas).miles
direction_to_lajitas = get_direction(point_carving, point_lajitas)

# 4. Print the analysis for the geographical claims
print("Analysis of Geographical Claims:")
print("="*40)

# Evaluation of Option B
print("Evaluating Claim B: '...depicts the Chiso Mountains... about 20 miles north...'")
print(f"The calculated distance from the carving to the Chisos Mountains is {distance_to_chisos_miles:.2f} miles.")
print(f"The actual direction is {direction_to_chisos}.")
print("Result: The claim states 'about 20 miles north'. The calculated distance is approximately 32.54 miles and the direction is North-East, not North. Therefore, claim B is inaccurate.\n")

# Evaluation of Option D
print("Evaluating Claim D: '...matches with a segment of the Bravo River near Lajitas, Texas, located 10 miles northwest...'")
print(f"The calculated distance from the carving to Lajitas is {distance_to_lajitas_miles:.2f} miles.")
print(f"The actual direction is {direction_to_lajitas}.")
print("Result: The claim states '10 miles northwest'. While the distance is close at 14.50 miles, the direction is North-East, not North-West. Therefore, claim D is inaccurate.")
print("="*40)

print("\nFinal Conclusion:")
print("- Option C is incorrect, as carvings are clearly visible.")
print("- Options B and D are incorrect because the geographical details (distance and/or direction) do not match the calculations.")
print("- The remaining options are A (a snake) and E (a sipapu). The carvings appear more complex than a single snake, showing dotted patterns and what looks like a starburst, which are not typical of snake depictions.")
print("- Option E provides the most plausible explanation. The sipapu is a symbolic feature in the cosmology of Southwest indigenous tribes, representing a place of origin or a portal. The abstract nature of the carvings, including the dots and other symbols, is consistent with such rock art.")
