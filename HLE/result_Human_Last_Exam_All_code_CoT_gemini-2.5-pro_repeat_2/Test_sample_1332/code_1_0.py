def calculate_minimum_storage():
    """
    Calculates the minimum storage space in bytes for a single record
    in the FlightRecords table.
    """
    # 1. Analyze flightNumber (integer, PRIMARY KEY)
    # A standard integer takes 4 bytes. As a PRIMARY KEY, it cannot be NULL.
    flightNumber_min_bytes = 4

    # 2. Analyze ticketCost (integer, nullable)
    # For minimum storage, we assume this field is NULL.
    # The storage for the data of a NULL value is 0 bytes.
    ticketCost_min_bytes = 0

    # 3. Analyze arrivalCity (text, nullable)
    # For minimum storage, we assume this field is NULL. Storage is 0 bytes.
    arrivalCity_min_bytes = 0

    # 4. Analyze departureCity (text, nullable)
    # For minimum storage, we assume this field is NULL. Storage is 0 bytes.
    departureCity_min_bytes = 0

    # 5. Calculate the total minimum size
    total_min_bytes = (
        flightNumber_min_bytes
        + ticketCost_min_bytes
        + arrivalCity_min_bytes
        + departureCity_min_bytes
    )

    print("Calculating the minimum storage space for one record:")
    print(f"- flightNumber (integer, NOT NULL): {flightNumber_min_bytes} bytes")
    print(f"- ticketCost (integer, NULLable): {ticketCost_min_bytes} bytes (when NULL)")
    print(f"- arrivalCity (text, NULLable): {arrivalCity_min_bytes} bytes (when NULL)")
    print(f"- departureCity (text, NULLable): {departureCity_min_bytes} bytes (when NULL)")
    print("\nFinal equation for minimum storage:")
    print(f"{flightNumber_min_bytes} + {ticketCost_min_bytes} + {arrivalCity_min_bytes} + {departureCity_min_bytes} = {total_min_bytes}")
    print(f"\nThe minimum storage space required is {total_min_bytes} bytes.")
    print("<<<4>>>")

calculate_minimum_storage()