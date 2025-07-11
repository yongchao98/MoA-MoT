def calculate_min_record_size():
    """
    Calculates the minimum storage space in bytes for a single record
    in the FlightRecords table.
    """
    # Size of a standard integer in bytes. This field cannot be NULL.
    flight_number_size = 4
    
    # The other fields can be NULL. For minimum size, we assume they are all NULL.
    # The data storage space for a NULL value is 0.
    ticket_cost_size_min = 0
    arrival_city_size_min = 0
    departure_city_size_min = 0
    
    # The total minimum size is the sum of the minimum sizes of each field.
    total_min_size = (
        flight_number_size +
        ticket_cost_size_min +
        arrival_city_size_min +
        departure_city_size_min
    )
    
    print("Minimum Storage Calculation:")
    print(f"flightNumber (integer, NOT NULL): {flight_number_size} bytes")
    print(f"ticketCost (nullable integer, set to NULL): {ticket_cost_size_min} bytes")
    print(f"arrivalCity (nullable text, set to NULL): {arrival_city_size_min} bytes")
    print(f"departureCity (nullable text, set to NULL): {departure_city_size_min} bytes")
    print("-" * 25)
    print(f"Total Minimum Size = {flight_number_size} + {ticket_cost_size_min} + {arrival_city_size_min} + {departure_city_size_min} = {total_min_size} bytes")

calculate_min_record_size()