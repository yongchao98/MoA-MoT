def calculate_min_record_size():
    """
    Calculates the minimum storage space in bytes for a single record
    in the FlightRecords table.

    Assumptions:
    - A standard integer takes 4 bytes.
    - NULL values for data take 0 bytes of storage for the data itself.
    - A null bitmap is used to track which nullable fields are NULL.
      Each nullable field requires one bit in the bitmap.
    - The calculation does not include database-specific overhead like
      record headers, which are constant for every row.
    """

    # 1. Storage for non-nullable fields
    # flightNumber is an integer PRIMARY KEY, so it cannot be NULL.
    size_flight_number = 4  # bytes
    print(f"Minimum size of non-nullable field (flightNumber as integer): {size_flight_number} bytes")

    # 2. Storage for null-tracking mechanism
    # There are 3 nullable fields: ticketCost, arrivalCity, departureCity.
    num_nullable_fields = 3
    # We need 1 bit per nullable field to track if it's NULL.
    # Total bits needed = 3.
    # This must be stored in a whole number of bytes. 3 bits fit in 1 byte.
    size_null_bitmap = 1  # bytes
    print(f"Minimum size for null tracking ({num_nullable_fields} fields require a 3-bit map, stored in 1 byte): {size_null_bitmap} byte")

    # 3. Calculate total minimum size
    # In the minimum-sized record, all nullable fields are NULL and their data size is 0.
    # The total size is the sum of the non-nullable data and the null bitmap.
    total_min_size = size_flight_number + size_null_bitmap

    print("\nFinal Calculation:")
    print(f"Total Minimum Storage = (size of non-nullable data) + (size of null bitmap)")
    print(f"Total Minimum Storage = {size_flight_number} + {size_null_bitmap} = {total_min_size} bytes")


if __name__ == "__main__":
    calculate_min_record_size()
