def calculate_min_record_size():
    """
    Calculates the minimum storage space in bytes for a record in the
    FlightRecords table.

    Assumptions:
    - An 'integer' data type occupies 4 bytes.
    - NULL values do not consume data space for their field.
    - A null bitmap is used to track NULL fields, requiring 1 bit per
      nullable column, with the total rounded up to the nearest byte.
    """

    # Size of the non-nullable PRIMARY KEY field 'flightNumber'
    size_flight_number = 4  # A standard integer is 4 bytes.

    # Number of fields that can be NULL
    nullable_fields_count = 3  # (ticketCost, arrivalCity, departureCity)

    # Calculate the size of the null bitmap overhead.
    # We need 1 bit for each of the 3 nullable fields.
    # Since storage is byte-aligned, 3 bits require 1 byte.
    # In Python, we can calculate this as ceil(bits / 8).
    import math
    size_null_bitmap_overhead = math.ceil(nullable_fields_count / 8.0)

    # The data size for nullable fields is 0 when they are NULL
    size_of_null_data = 0

    # The total minimum size is the sum of the non-nullable fields' data
    # plus the null bitmap overhead.
    total_min_size = size_flight_number + size_null_bitmap_overhead

    print("--- Minimum Record Storage Calculation ---")
    print(f"1. Size of 'flightNumber' (integer, non-nullable): {size_flight_number} bytes")
    print(f"2. Data size of the 3 nullable fields (when NULL): {size_of_null_data} bytes")
    print(f"3. Overhead for null bitmap ({nullable_fields_count} fields): {size_null_bitmap_overhead} byte")
    print("------------------------------------------")
    print("Total minimum size = (Size of non-nullable fields) + (Null overhead)")
    # The final equation with each number explicitly printed
    print(f"Final Equation: {size_flight_number} + {size_null_bitmap_overhead} = {int(total_min_size)}")
    print("------------------------------------------")

calculate_min_record_size()