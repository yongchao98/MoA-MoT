def calculate_min_record_size():
    """
    Calculates the minimum storage space in bytes for the data payload of a single
    record in the specified table structure.
    """
    
    # Standard size for an integer data type in most database systems.
    size_integer = 4
    
    # The size contribution of a NULL value to the data payload is zero.
    size_when_null = 0
    
    # 1. flightNumber (integer, PRIMARY KEY)
    # As a PRIMARY KEY, this field cannot be NULL. It will always take up the space
    # required for an integer.
    min_size_flightNumber = size_integer
    
    # 2. ticketCost (integer, nullable)
    # To find the minimum record size, we assume this field is NULL.
    min_size_ticketCost = size_when_null
    
    # 3. arrivalCity (text, nullable)
    # To find the minimum record size, we assume this field is NULL.
    min_size_arrivalCity = size_when_null
    
    # 4. departureCity (text, nullable)
    # To find the minimum record size, we assume this field is NULL.
    min_size_departureCity = size_when_null
    
    # Calculate the total minimum size by summing the minimum size of each field.
    total_min_data_size = (min_size_flightNumber +
                             min_size_ticketCost +
                             min_size_arrivalCity +
                             min_size_departureCity)
                             
    print("Calculating the minimum storage space for the record's data payload:")
    print("-" * 60)
    print(f"Minimum size for 'flightNumber' (integer, NOT NULL): {min_size_flightNumber} bytes")
    print(f"Minimum size for 'ticketCost' (integer, when NULL):  {min_size_ticketCost} bytes")
    print(f"Minimum size for 'arrivalCity' (text, when NULL):   {min_size_arrivalCity} bytes")
    print(f"Minimum size for 'departureCity' (text, when NULL):  {min_size_departureCity} bytes")
    print("-" * 60)
    
    # Print the final equation as requested
    print(f"Final Equation: {min_size_flightNumber} + {min_size_ticketCost} + {min_size_arrivalCity} + {min_size_departureCity} = {total_min_data_size}")
    
    print("\nNote: This is the size of the data payload. Actual storage on disk includes additional system-specific overhead per record.")

calculate_min_record_size()