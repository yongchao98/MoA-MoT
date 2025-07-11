def calculate_minimum_record_size():
    """
    Calculates the minimum storage space in bytes for a single record
    in the specified table structure.

    The calculation is based on these assumptions:
    1. A fixed-size record header is required for metadata (e.g., transaction info, null bitmap).
       A typical size for this header is 24 bytes.
    2. 'flightNumber' is an integer PRIMARY KEY, so it cannot be null and requires 4 bytes.
    3. For minimum size, all nullable fields ('ticketCost', 'arrivalCity', 'departureCity')
       are assumed to be NULL. A NULL value occupies 0 bytes in the data area, as its
       status is tracked in the record header.
    """
    # Size in bytes for the record's metadata header
    record_header_size = 24

    # Size in bytes for each field under minimum storage conditions
    flight_number_size = 4  # integer, NOT NULL
    ticket_cost_size = 0      # integer, nullable, assumed NULL
    arrival_city_size = 0     # text, nullable, assumed NULL
    departure_city_size = 0   # text, nullable, assumed NULL

    # Calculate the total minimum size
    total_size = (record_header_size +
                  flight_number_size +
                  ticket_cost_size +
                  arrival_city_size +
                  departure_city_size)

    # Print the breakdown of the calculation
    print("Minimum record storage space calculation:")
    print(f"Record Header: {record_header_size} bytes")
    print(f"flightNumber (integer, non-null): {flight_number_size} bytes")
    print(f"ticketCost (nullable, assumed null): {ticket_cost_size} bytes")
    print(f"arrivalCity (nullable, assumed null): {arrival_city_size} bytes")
    print(f"departureCity (nullable, assumed null): {departure_city_size} bytes")
    print("---")
    print(f"Total = {record_header_size} + {flight_number_size} + {ticket_cost_size} + {arrival_city_size} + {departure_city_size} = {total_size} bytes")

calculate_minimum_record_size()