#
# This script calculates the minimum storage space in bytes for a single record
# based on the SQLite record format.
#

# --- Step 1: Calculate the size of the record header ---

# To minimize size, we choose values that have the smallest representation.
# - flightNumber: value 0 (special type code, 1 byte in header, 0 bytes in body)
# - Other fields: NULL (special type code, 1 byte in header, 0 bytes in body)

# The size of the type code for each column in the header (in bytes)
flight_number_type_size = 1
ticket_cost_type_size = 1
arrival_city_type_size = 1
departure_city_type_size = 1

# The header also contains a varint specifying the total size of the type codes.
# The total size is 4 bytes (1+1+1+1), which can be stored in a 1-byte varint.
header_size_varint_size = 1

# The total header size is the sum of the size-prefix varint and all the type codes.
total_header_size = header_size_varint_size + flight_number_type_size + ticket_cost_type_size + arrival_city_type_size + departure_city_type_size

# --- Step 2: Calculate the size of the record body ---

# For the chosen minimal values, the data stored in the body is:
flight_number_data_size = 0  # 0 bytes for the integer value 0
ticket_cost_data_size = 0    # 0 bytes for a NULL value
arrival_city_data_size = 0   # 0 bytes for a NULL value
departure_city_data_size = 0 # 0 bytes for a NULL value

# The total body size is the sum of the data sizes.
total_body_size = flight_number_data_size + ticket_cost_data_size + arrival_city_data_size + departure_city_data_size

# --- Step 3: Calculate the total minimum record size ---

total_record_size = total_header_size + total_body_size

print("Calculating the total minimum storage space for a record:")
print("---------------------------------------------------------")
print(f"Total Header Size = (Size Varint) + (Type Codes)")
print(f"Total Header Size = {header_size_varint_size} + ({flight_number_type_size} + {ticket_cost_type_size} + {arrival_city_type_size} + {departure_city_type_size}) = {total_header_size} bytes")
print("")
print(f"Total Body Size = (flightNumber) + (ticketCost) + (arrivalCity) + (departureCity)")
print(f"Total Body Size = {flight_number_data_size} + {ticket_cost_data_size} + {arrival_city_data_size} + {departure_city_data_size} = {total_body_size} bytes")
print("---------------------------------------------------------")
print("Final Calculation:")
print(f"Total Record Size = Total Header Size + Total Body Size")
print(f"Total Record Size = {total_header_size} + {total_body_size} = {total_record_size} bytes")
<<<5>>>