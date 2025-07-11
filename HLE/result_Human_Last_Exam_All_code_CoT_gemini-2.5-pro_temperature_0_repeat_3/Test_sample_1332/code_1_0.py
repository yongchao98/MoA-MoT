# Plan:
# 1. Define the size in bytes for each field in its minimal state.
#    - A non-nullable integer ('flightNumber') always takes 4 bytes.
#    - A nullable field ('ticketCost', 'arrivalCity', 'departureCity')
#      takes 0 bytes for data storage when it is NULL.
# 2. Sum the sizes to find the total minimum storage space for the record's data.
# 3. Print the calculation and the final result.

# Size of the non-nullable integer primary key
flightNumber_size = 4

# Size of nullable fields when they are NULL is 0
ticketCost_min_size = 0
arrivalCity_min_size = 0
departureCity_min_size = 0

# Calculate the total minimum size
total_min_size = flightNumber_size + ticketCost_min_size + arrivalCity_min_size + departureCity_min_size

# Print the final equation and result
print(f"The minimum storage space is calculated by summing the sizes of each field in its smallest state.")
print(f"Minimum size = size of flightNumber + size of ticketCost (NULL) + size of arrivalCity (NULL) + size of departureCity (NULL)")
print(f"Minimum size = {flightNumber_size} bytes + {ticketCost_min_size} bytes + {arrivalCity_min_size} bytes + {departureCity_min_size} bytes = {total_min_size} bytes")