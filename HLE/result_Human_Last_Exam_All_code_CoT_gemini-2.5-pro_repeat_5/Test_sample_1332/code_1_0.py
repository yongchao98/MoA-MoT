# plan:
# 1. Define the storage size in bytes for the non-nullable field.
#    - A standard integer takes 4 bytes.
flightNumber_size = 4

# 2. Define the storage size for the nullable fields when they are NULL.
#    - When a field is NULL, its data storage cost is 0 bytes.
ticketCost_size = 0
arrivalCity_size = 0
departureCity_size = 0

# 3. Calculate the total minimum size by summing the individual sizes.
total_min_size = flightNumber_size + ticketCost_size + arrivalCity_size + departureCity_size

# 4. Print the calculation and the final result.
print("Minimum storage space calculation:")
print(f"flightNumber (integer, NOT NULL): {flightNumber_size} bytes")
print(f"ticketCost (integer, NULL): {ticketCost_size} bytes")
print(f"arrivalCity (text, NULL): {arrivalCity_size} bytes")
print(f"departureCity (text, NULL): {departureCity_size} bytes")
print(f"\nTotal minimum size = {flightNumber_size} + {ticketCost_size} + {arrivalCity_size} + {departureCity_size} = {total_min_size} bytes")
