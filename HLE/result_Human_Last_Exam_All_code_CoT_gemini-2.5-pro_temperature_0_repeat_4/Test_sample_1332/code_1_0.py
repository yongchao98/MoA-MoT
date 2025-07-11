def calculate_min_record_size():
    """
    Calculates the minimum storage space in bytes for a single record
    in the specified table structure.
    """

    # --- Step 1: Define component sizes in bytes ---
    # Based on common database systems (e.g., PostgreSQL)

    # Size of a standard integer.
    integer_size = 4

    # Size of a field when it is NULL. It consumes no data space,
    # as its null status is tracked in the record header's null bitmap.
    size_when_null = 0

    # Fixed overhead per record for metadata (header, transaction info, null bitmap, etc.).
    # A typical value is 24 bytes for robust database systems.
    record_overhead = 24

    # --- Step 2: Determine the size of each field for a minimal record ---

    # flightNumber (integer, PRIMARY KEY): This field cannot be NULL.
    flight_number_size = integer_size

    # ticketCost (integer, nullable): To minimize size, we assume it's NULL.
    ticket_cost_size = size_when_null

    # arrivalCity (text, nullable): Assumed to be NULL for minimum size.
    arrival_city_size = size_when_null

    # departureCity (text, nullable): Assumed to be NULL for minimum size.
    departure_city_size = size_when_null

    # --- Step 3: Calculate the total minimum size ---
    total_size = (record_overhead +
                  flight_number_size +
                  ticket_cost_size +
                  arrival_city_size +
                  departure_city_size)

    # --- Step 4: Print the explanation and final equation ---
    print("Calculating the minimum storage space for one record:")
    print(f"- Record Overhead (metadata, header): {record_overhead} bytes")
    print(f"- flightNumber (integer, non-nullable): {flight_number_size} bytes")
    print(f"- ticketCost (integer, nullable): {ticket_cost_size} bytes (when NULL)")
    print(f"- arrivalCity (text, nullable): {arrival_city_size} bytes (when NULL)")
    print(f"- departureCity (text, nullable): {departure_city_size} bytes (when NULL)")
    print("\nFinal Equation (in bytes):")
    print(f"{record_overhead} + {flight_number_size} + {ticket_cost_size} + {arrival_city_size} + {departure_city_size} = {total_size}")

calculate_min_record_size()
<<<28>>>