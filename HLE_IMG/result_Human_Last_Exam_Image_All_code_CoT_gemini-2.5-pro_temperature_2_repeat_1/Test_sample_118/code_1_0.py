# This script represents the solution to the pathfinding problem.
# The optimal path's details have been determined by analyzing the Zurich Tram/Bus map
# according to the specified criteria.

# --- Path Components ---

# S stands for the start station, Billoweg.
start_symbol = "S"

# The first leg of the journey is on line 7.
line_1 = 7

# The exchange happens at Bahnhofstrasse/HB.
exchange_station_1 = "Bahnhofstrasse/HB"

# The second leg of the journey is on line 14.
line_2 = 14

# E stands for the end station, Laubiweg.
end_symbol = "E"

# The total number of stations on this path, including start and end.
total_stations = 19

# --- Format and Print the Result ---

# The final answer is formatted as specified in the problem description.
# The path describes going from S to E via the exchange station, using the given lines.
# The number after the semicolon indicates the total station count.
final_path_description = f"{start_symbol} - {line_1} - {exchange_station_1} - {line_2} - {end_symbol}; {total_stations}"

print(final_path_description)