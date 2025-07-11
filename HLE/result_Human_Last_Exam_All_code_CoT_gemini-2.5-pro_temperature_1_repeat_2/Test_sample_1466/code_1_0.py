import geopy.distance

def find_fastest_route():
    """
    Calculates the estimated walking time for five different routes from
    Guildhall to St Paul's Cathedral, accounting for a road closure.
    """
    # Step 1 & 2: Define coordinates for key locations and waypoints.
    locations = {
        'Guildhall': (51.5155, -0.0918),
        'St Pauls Cathedral': (51.5138, -0.0984),
        'A_Gresham_Foster': (51.5155, -0.0950),     # Route A: Intersection of Gresham St and Foster Ln
        'A_Cheapside_NewChange': (51.5142, -0.0975),# Route A: Intersection of Cheapside and New Change
        'B_E_MuseumOfLondon': (51.5173, -0.0963),   # Route B/E: Represents Aldersgate St roundabout
        'B_E_Cheapside_StMartins': (51.5144, -0.0984),# Route B/E: Intersection of Cheapside and St Martin's Le Grand
        'C_BankOfEngland': (51.5136, -0.0886),      # Route C: Represents Lothbury/Princes St area
        'C_MansionHouseStation': (51.5123, -0.0924),# Route C: On Queen Victoria St
        'D_CannonStStation': (51.5113, -0.0905),    # Route D: A key point on the southern detour
        'D_MillenniumBridge': (51.5106, -0.0976),   # Route D: Another key point on the southern detour
        'E_WoodSt_Barbican': (51.5175, -0.0934),    # Route E: Northernmost point of the detour
    }

    routes = {
        'A': ['Guildhall', 'A_Gresham_Foster', 'A_Cheapside_NewChange', 'St Pauls Cathedral'],
        'B': ['Guildhall', 'B_E_MuseumOfLondon', 'B_E_Cheapside_StMartins', 'St Pauls Cathedral'],
        'C': ['Guildhall', 'C_BankOfEngland', 'C_MansionHouseStation', 'St Pauls Cathedral'],
        'D': ['Guildhall', 'D_CannonStStation', 'D_MillenniumBridge', 'St Pauls Cathedral'],
        'E': ['Guildhall', 'E_WoodSt_Barbican', 'B_E_MuseumOfLondon', 'B_E_Cheapside_StMartins', 'St Pauls Cathedral']
    }

    # Constants for calculation
    WALKING_SPEED_KMH = 5.0
    CIRCUITY_FACTOR = 1.4  # Factor to estimate road distance from straight-line distance

    route_times = {}

    print("Calculating the fastest walking route from Guildhall to St Paul's Cathedral...")
    print(f"Assuming an average walking speed of {WALKING_SPEED_KMH} km/h and a city circuity factor of {CIRCUITY_FACTOR}.\n")

    # Step 3, 4, 5: Calculate distance and time for each route
    for name, waypoints in routes.items():
        total_dist_km = 0
        dist_segments = []
        
        waypoint_coords = [locations[wp] for wp in waypoints]
        
        for i in range(len(waypoint_coords) - 1):
            dist = geopy.distance.great_circle(waypoint_coords[i], waypoint_coords[i+1]).km
            total_dist_km += dist
            dist_segments.append(f"{dist:.3f}")

        road_dist_km = total_dist_km * CIRCUITY_FACTOR
        time_min = (road_dist_km / WALKING_SPEED_KMH) * 60
        route_times[name] = time_min

        # Output the calculation for each route
        dist_sum_str = " + ".join(dist_segments)
        print(f"Route {name}:")
        print(f"  Straight-line distance (km): {dist_sum_str} = {total_dist_km:.3f} km")
        print(f"  Estimated road distance (km): {total_dist_km:.3f} * {CIRCUITY_FACTOR} = {road_dist_km:.3f} km")
        print(f"  Estimated time (minutes): {road_dist_km:.3f} km / {WALKING_SPEED_KMH} km/h * 60 = {time_min:.2f} minutes\n")

    # Step 6: Find and announce the fastest route
    fastest_route = min(route_times, key=route_times.get)
    fastest_time = route_times[fastest_route]

    print(f"--- Conclusion ---")
    print(f"Based on the calculations, the fastest option is Route {fastest_route} with an estimated time of {fastest_time:.2f} minutes.")

if __name__ == '__main__':
    find_fastest_route()