# Path details based on the analysis
start_station_name = "Billoweg"
end_station_name = "Laubiweg"
exchange_station = "Milchbuck"
first_line = 7
second_line = 9
total_stations = 17

# Format the output string as requested.
# The 'S' and 'E' are placeholders for the start and end stations.
path_description = f"S - {first_line} - {exchange_station} - {second_line} - E; {total_stations}"

print(path_description)

# We are asked to output each number in the final equation.
# The numbers are the line numbers and the total station count.
print(f"First line: {first_line}")
print(f"Second line: {second_line}")
print(f"Total stations: {total_stations}")