import functools

class Path:
    """A class to represent a path with multiple legs."""
    def __init__(self, legs):
        self.legs = legs
        self.exchanges = len(legs) - 1
        self.station_count = self._calculate_station_count()
        self.exchange_points = self._get_exchange_points()
        self.total_stations_at_first_exchange = self._get_stations_at_first_exchange()

    def _calculate_station_count(self):
        if not self.legs:
            return 0
        total_stations = len(self.legs[0])
        for i in range(1, len(self.legs)):
            # Add stations from the new leg, excluding the transfer station
            total_stations += len(self.legs[i]) - 1
        return total_stations
    
    def _get_exchange_points(self):
        points = []
        for i in range(len(self.legs) - 1):
            points.append(self.legs[i][-1])
        return points
        
    def _get_stations_at_first_exchange(self):
        if not self.legs or len(self.legs) < 2:
            return float('inf')
        return len(self.legs[0])

    def __repr__(self):
        path_str = "S"
        for i, leg in enumerate(self.legs):
            line, stations = leg.line, leg.stations
            if i > 0:
                 path_str += f" - {leg.stations[0]}"
            path_str += f" - {line}"
        path_str += f" - E; {self.station_count}"
        return path_str

class Leg:
    """A class to represent one segment of a journey on a single line."""
    def __init__(self, line, stations):
        self.line = line
        self.stations = stations
    
    def __len__(self):
        return len(self.stations)
    
    def __getitem__(self, key):
        return self.stations[key]


def compare_paths(path1, path2):
    """Comparison function to sort paths based on the given criteria."""
    # A) Fewer exchanges
    if path1.exchanges != path2.exchanges:
        return path1.exchanges - path2.exchanges
    
    # B) Fewer stations
    if path1.station_count != path2.station_count:
        return path1.station_count - path2.station_count
    
    # C) Later exchanges (modified interpretation)
    # We compare the number of stations traveled before the first exchange.
    if path1.total_stations_at_first_exchange != path2.total_stations_at_first_exchange:
        # We want the larger value, so we reverse the subtraction order.
        return path2.total_stations_at_first_exchange - path1.total_stations_at_first_exchange
        
    # Final tie-breaker: lower line number on the second leg
    if len(path1.legs) > 1 and len(path2.legs) > 1:
        if path1.legs[1].line != path2.legs[1].line:
            return path1.legs[1].line - path2.legs[1].line
            
    return 0

def solve_zurich_path():
    """
    Finds the optimal path from Billoweg to Laubiweg based on predefined rules.
    """
    # Define stations for each potential path segment
    leg_72_to_paradeplatz = Leg(72, ['Billoweg', 'Brunaustr.', 'Hügelstr.', 'Bhf Enge', 'Tunnelstr.', 'Stockerstr.', 'Paradeplatz'])
    leg_11_from_paradeplatz = Leg(11, ['Paradeplatz', 'Bahnhofstrasse/HB', 'Bahnhofquai/HB', 'Stampfenbachplatz', 'Beckenhof', 'Schaffhauserplatz', 'Guggachstr.', 'Laubiweg'])

    leg_72_to_stampfenbachplatz = Leg(72, leg_72_to_paradeplatz.stations + ['Rennweg', 'Bahnhofplatz/HB', 'Stampfenbachplatz'])
    leg_11_15_from_stampfenbachplatz = ['Stampfenbachplatz', 'Beckenhof', 'Schaffhauserplatz', 'Guggachstr.', 'Laubiweg']
    
    leg_72_to_beckenhof = Leg(72, leg_72_to_stampfenbachplatz.stations + ['Beckenhof'])
    leg_11_15_from_beckenhof = ['Beckenhof', 'Schaffhauserplatz', 'Guggachstr.', 'Laubiweg']

    leg_72_to_schaffhauserplatz = Leg(72, leg_72_to_beckenhof.stations + ['Schaffhauserplatz'])
    leg_11_15_from_schaffhauserplatz = ['Schaffhauserplatz', 'Guggachstr.', 'Laubiweg']

    # Create all 1-exchange candidate paths
    candidate_paths = [
        Path([leg_72_to_paradeplatz, leg_11_from_paradeplatz]),
        
        Path([leg_72_to_stampfenbachplatz, Leg(11, leg_11_15_from_stampfenbachplatz)]),
        Path([leg_72_to_stampfenbachplatz, Leg(15, leg_11_15_from_stampfenbachplatz)]),
        
        Path([leg_72_to_beckenhof, Leg(11, leg_11_15_from_beckenhof)]),
        Path([leg_72_to_beckenhof, Leg(15, leg_11_15_from_beckenhof)]),

        Path([leg_72_to_schaffhauserplatz, Leg(11, leg_11_15_from_schaffhauserplatz)]),
        Path([leg_72_to_schaffhauserplatz, Leg(15, leg_11_15_from_schaffhauserplatz)]),
    ]

    # Sort the paths using the defined criteria
    candidate_paths.sort(key=functools.cmp_to_key(compare_paths))
    
    # The best path is the first one in the sorted list
    best_path = candidate_paths[0]

    # Format the final answer string
    exchange_station = best_path.exchange_points[0]
    first_line = best_path.legs[0].line
    second_line = best_path.legs[1].line
    station_count = best_path.station_count

    final_answer = f"S - {first_line} - {exchange_station} - {second_line} - E; {station_count}"
    print(final_answer)

solve_zurich_path()