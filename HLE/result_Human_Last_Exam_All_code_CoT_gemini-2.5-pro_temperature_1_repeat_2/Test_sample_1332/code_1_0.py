def calculate_min_record_size():
    """
    Calculates the minimum storage space in bytes for a record in the FlightRecords table.
    """
    # Standard size in bytes for an integer data type.
    size_integer = 4

    # For minimum record size, nullable fields are assumed to be NULL,
    # which contributes 0 bytes to the data storage.
    size_null_field = 0

    # 1. flightNumber: As a PRIMARY KEY, it cannot be NULL.
    min_size_flightNumber = size_integer

    # 2. ticketCost: Can be NULL.
    min_size_ticketCost = size_null_field

    # 3. arrivalCity: Can be NULL.
    min_size_arrivalCity = size_null_field

    # 4. departureCity: Can be NULL.
    min_size_departureCity = size_null_field

    # Total minimum size is the sum of the minimum size of each field.
    total_min_size = (min_size_flightNumber +
                      min_size_ticketCost +
                      min_size_arrivalCity +
                      min_size_departureCity)

    print("Minimum storage space calculation:")
    print(f"flightNumber (integer, NOT NULL): {min_size_flightNumber} bytes")
    print(f"ticketCost (integer, NULL): {min_size_ticketCost} bytes")
    print(f"arrivalCity (text, NULL): {min_size_arrivalCity} bytes")
    print(f"departureCity (text, NULL): {min_size_departureCity} bytes")
    print("-" * 30)
    print(f"Total minimum size = {min_size_flightNumber} + {min_size_ticketCost} + {min_size_arrivalCity} + {min_size_departureCity} = {total_min_size} bytes")

calculate_min_record_size()
<<<4>>>