def calculate_exposure_time():
    """
    This function calculates the duration the PDA media was exposed to room air.
    """
    # The media was exposed from 7 am to 1 pm.
    # In 24-hour format, this is from 7:00 to 13:00.
    start_time_24hr = 7
    discovery_time_24hr = 13

    # Calculate the duration
    duration = discovery_time_24hr - start_time_24hr

    # Print the equation and result
    print("The laboratory discovered the media had been exposed to room air for a significant duration.")
    print(f"The calculation for the exposure time is: {discovery_time_24hr} - {start_time_24hr} = {duration} hours.")

calculate_exposure_time()