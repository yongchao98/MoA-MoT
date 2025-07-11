def calculate_minimum_storage():
    """
    Calculates the minimum storage space in bytes for a record in the specified table.
    """
    
    # 1. Size of the non-nullable primary key field
    # 'flightNumber' is an integer, which typically requires 4 bytes.
    # As a PRIMARY KEY, it cannot be NULL.
    flight_number_size = 4
    
    # 2. Size of nullable fields at their minimum
    # 'ticketCost', 'arrivalCity', and 'departureCity' can all be NULL.
    # When a field is NULL, its data storage cost is 0 bytes.
    # The information about them being NULL is stored separately.
    nullable_fields_data_size = 0
    
    # 3. Size of the metadata to track NULLs (Null Bitmap)
    # There are 3 nullable fields. We need 1 bit for each to flag if it's NULL.
    # This requires 3 bits. Since the smallest unit of storage is a byte (8 bits),
    # we need 1 byte to store this information.
    null_bitmap_size = 1
    
    # 4. Calculate the total minimum size
    total_min_size = flight_number_size + nullable_fields_data_size + null_bitmap_size
    
    print("Calculating the minimum storage space for a single record:")
    print("-" * 50)
    print("1. 'flightNumber' (integer, PRIMARY KEY): This field cannot be null.")
    print(f"   - Storage for integer = {flight_number_size} bytes")
    print("\n2. 'ticketCost', 'arrivalCity', 'departureCity' (nullable fields):")
    print("   - For minimum record size, these are assumed to be NULL.")
    print(f"   - Storage for NULL data = {nullable_fields_data_size} bytes")
    print("\n3. Null Bitmap: The database must track which of the 3 fields are NULL.")
    print("   - This requires 3 bits, which takes up 1 byte of storage.")
    print(f"   - Storage for null bitmap = {null_bitmap_size} byte")
    print("-" * 50)
    print("Total Minimum Storage = (size of flightNumber) + (size of null bitmap)")
    print(f"Equation: {flight_number_size} + {null_bitmap_size} = {total_min_size} bytes")

calculate_minimum_storage()
<<<5>>>