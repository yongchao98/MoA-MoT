import collections

def calculate_periods():
    """
    Calculates the 'period' of tiles defined by lists of points.
    The period is interpreted as the number of points defining the tile.
    """
    point_sets = collections.OrderedDict()
    point_sets["Tile 1"] = [13, 31, 23]
    point_sets["Tile 2"] = [10, 4, 23, 31]
    point_sets["Tile 3"] = [5, 15, 17, 19, 21, 7]
    point_sets["Tile 4"] = [4, 5, 14, 23, 18, 19, 21, 22, 31, 30, 9, 8, 13]
    
    periods = []
    
    # Calculate and show the equation for each period
    for name, points in point_sets.items():
        period = len(points)
        periods.append(period)
        print(f"Period of {name} ({points}) = {period}")
        
    print("\nThe four periods separated by ',' are:")
    # Print the final result in the requested format
    print(','.join(map(str, periods)))

calculate_periods()