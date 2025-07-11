def calculate_min_record_size():
    """
    Calculates the minimum storage space in bytes for a single record
    in the FlightRecords table.
    """
    # Size of an integer type in bytes. A primary key cannot be null.
    flight_number_size = 4

    # The other fields can be null, so their minimum data storage is 0 bytes.
    ticket_cost_min_size = 0
    arrival_city_min_size = 0
    departure_city_min_size = 0

    # The total minimum size is the sum of the minimum sizes of each field.
    total_min_size = (flight_number_size + ticket_cost_min_size +
                      arrival_city_min_size + departure_city_min_size)

    print("To find the minimum storage space, we sum the minimum size of each field:")
    print(f"flightNumber (integer, not null): {flight_number_size} bytes")
    print(f"ticketCost (integer, nullable): {ticket_cost_min_size} bytes (when NULL)")
    print(f"arrivalCity (text, nullable): {arrival_city_min_size} bytes (when NULL)")
    print(f"departureCity (text, nullable): {departure_city_min_size} bytes (when NULL)")
    print("-" * 30)
    print(f"Equation: {flight_number_size} + {ticket_cost_min_size} + {arrival_city_min_size} + {departure_city_min_size} = {total_min_size}")
    print(f"The minimum storage space required for the record's data is {total_min_size} bytes.")

calculate_min_record_size()