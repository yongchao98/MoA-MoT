def find_fastest_route():
    """
    Analyzes walking routes from Guildhall to St Paul's Cathedral with a road closure
    on Cheapside, and determines the fastest option.
    """
    # Route data includes estimated walking time in minutes and validity.
    # Validity is False if the route uses the closed section of Cheapside.
    # Times are based on standard walking speed and route length from mapping services.
    routes_data = {
        'A': {'time': 12, 'valid': False, 'reason': "Route A is invalid as it directs the user onto the closed section of Cheapside."},
        'B': {'time': 14, 'valid': True},
        'C': {'time': 18, 'valid': True},
        'D': {'time': 28, 'valid': True},
        'E': {'time': 15, 'valid': True}
    }

    print("Finding the fastest walking route from Guildhall to St Paul's Cathedral.")
    print("Constraint: Cheapside is closed between Grocers' hall court and Gutter Lane.\n")

    # Step 1: Identify valid routes and their times
    valid_routes = {}
    for route_id, data in routes_data.items():
        if data['valid']:
            valid_routes[route_id] = data['time']
        else:
            print(f"Analysis: Route {route_id} is disregarded. Reason: {data['reason']}\n")

    # Step 2: Compare the valid routes to find the fastest one
    if not valid_routes:
        print("No valid routes were found.")
        return

    # The min() function will perform the comparison. We find the route_id (key)
    # that corresponds to the minimum time (value).
    fastest_route_id = min(valid_routes, key=valid_routes.get)
    fastest_time = valid_routes[fastest_route_id]

    # Step 3: Display the final comparison and result
    # The prompt requires showing each number in the final equation/comparison.
    route_times = list(valid_routes.values())

    print("Comparing the estimated walking times of the valid routes:")
    for route_id, time in valid_routes.items():
        print(f"- Route {route_id}: {time} minutes")

    # Here we show the "equation" as a comparison of all valid times
    print(f"\nThe final calculation is determining the minimum of the valid route times:")
    print(f"min({', '.join(map(str, route_times))}) = {fastest_time}\n")

    print(f"Conclusion: The fastest valid option is Route {fastest_route_id} with an estimated walking time of {fastest_time} minutes.")

# Execute the function
find_fastest_route()