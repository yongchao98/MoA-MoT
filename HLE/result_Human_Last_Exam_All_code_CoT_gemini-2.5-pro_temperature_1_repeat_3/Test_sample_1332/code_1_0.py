def calculate_minimum_record_size():
    """
    Calculates the minimum storage space in bytes for a record in the FlightRecords table.
    This calculation is based on the size of the data types, assuming nullable fields are NULL.
    """
    # Standard size in bytes for an integer data type.
    integer_size = 4

    # The size contribution of a NULL value is 0 data bytes.
    null_size = 0

    # 1. flightNumber (integer, PRIMARY KEY): Cannot be NULL.
    flight_number_min_size = integer_size

    # 2. ticketCost (integer, nullable): Smallest size is when it's NULL.
    ticket_cost_min_size = null_size

    # 3. arrivalCity (text, nullable): Smallest size is when it's NULL.
    arrival_city_min_size = null_size

    # 4. departureCity (text, nullable): Smallest size is when it's NULL.
    departure_city_min_size = null_size

    # Sum the sizes for the total minimum record size.
    total_min_size = (flight_number_min_size + 
                      ticket_cost_min_size + 
                      arrival_city_min_size + 
                      departure_city_min_size)

    print("Calculating the minimum storage space for a single record:")
    print("This is achieved by setting all nullable fields to NULL.")
    print("-" * 50)
    print(f"Size of 'flightNumber' (integer, NOT NULL): {flight_number_min_size} bytes")
    print(f"Size of 'ticketCost' (nullable): {ticket_cost_min_size} bytes")
    print(f"Size of 'arrivalCity' (nullable): {arrival_city_min_size} bytes")
    print(f"Size of 'departureCity' (nullable): {departure_city_min_size} bytes")
    print("-" * 50)
    print("Final Equation:")
    print(f"{flight_number_min_size} + {ticket_cost_min_size} + {arrival_city_min_size} + {departure_city_min_size} = {total_min_size}")
    print(f"\nThe minimum storage space required is {total_min_size} bytes.")

calculate_minimum_record_size()