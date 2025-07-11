import math

def calculate_min_record_size():
    """
    Calculates the minimum storage space in bytes for a record in the specified table.
    This assumes nullable fields are set to NULL for the smallest possible size.
    """
    # 1. Size of non-nullable fields.
    # flightNumber is an integer PRIMARY KEY, so it's not nullable.
    # A standard integer takes 4 bytes.
    size_of_non_nullable_data = 4
    field_name = "flightNumber"
    print(f"Storage for the non-nullable field '{field_name}' (integer) = {size_of_non_nullable_data} bytes.")

    # 2. Identify nullable fields to calculate the null bitmap size.
    # The fields ticketCost, arrivalCity, and departureCity can be null.
    num_nullable_fields = 3
    print(f"Number of nullable fields ('ticketCost', 'arrivalCity', 'departureCity') = {num_nullable_fields}.")

    # 3. Calculate the size of the null bitmap.
    # It requires 1 bit per nullable field, rounded up to the nearest byte.
    size_of_null_bitmap = math.ceil(num_nullable_fields / 8.0)
    print(f"Storage for null bitmap = ceil({num_nullable_fields} / 8) = {int(size_of_null_bitmap)} byte(s).")

    # 4. For minimum size, the data storage for nullable fields is 0 bytes.
    print(f"Storage for data in nullable fields (when NULL) = 0 bytes.")

    # 5. The total minimum size is the sum of non-nullable data and the null bitmap.
    total_min_size = size_of_non_nullable_data + size_of_null_bitmap

    print("\n--- Final Calculation ---")
    print(f"Minimum Record Size = (Size of '{field_name}') + (Size of Null Bitmap)")
    print(f"                    = {size_of_non_nullable_data} + {int(size_of_null_bitmap)}")
    print(f"                    = {int(total_min_size)}")

calculate_min_record_size()