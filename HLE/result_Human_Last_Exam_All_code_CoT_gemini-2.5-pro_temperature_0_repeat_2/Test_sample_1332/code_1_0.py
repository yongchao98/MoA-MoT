def calculate_min_record_size():
    """
    Calculates the minimum storage space in bytes for a record in the specified table.

    The calculation is based on these assumptions for minimum size:
    1. A fixed-size header for each record (tuple). A standard size is 24 bytes.
    2. The primary key 'flightNumber' (integer) is not nullable and requires 4 bytes.
    3. All other fields are nullable and are set to NULL for the minimum size,
       meaning their data storage size is 0 bytes.
    """

    # Size in bytes for the fixed record header (includes space for NULL bitmap)
    header_size = 24

    # Size in bytes for the non-nullable 'flightNumber' field (integer)
    flightnumber_size = 4

    # For minimum record size, all nullable fields are NULL, so their data size is 0.
    nullable_fields_size = 0

    # Total minimum size is the sum of the header and the required data fields.
    total_min_size = header_size + flightnumber_size + nullable_fields_size

    print("Calculating the minimum storage space for one record:")
    print(f"Record Header Size: {header_size} bytes")
    print(f"flightNumber (integer, NOT NULL) Size: {flightnumber_size} bytes")
    print(f"Nullable Fields (all NULL) Size: {nullable_fields_size} bytes")
    print("-" * 30)
    print("Final Equation (in bytes):")
    print(f"{total_min_size} = {header_size} + {flightnumber_size}")

calculate_min_record_size()
<<<28>>>