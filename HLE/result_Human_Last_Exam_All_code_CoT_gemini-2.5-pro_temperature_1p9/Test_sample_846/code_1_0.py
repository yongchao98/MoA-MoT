import datetime

def get_azobenzene_table_state():
    """
    Determines and prints the state of the azobenzene-shaped table based on the time of day.
    This simulates the molecule's photoisomerization.
    - Daytime (with UV sunlight): Bent 'cis' form.
    - Nighttime (no UV light): Straight 'trans' form.
    """
    # Define approximate sunrise and sunset hours for the simulation
    # Sunrise is at 6 AM
    sunrise_hour = 6
    # Sunset is at 8 PM (20:00 in 24-hour format)
    sunset_hour = 20
    
    current_hour = datetime.datetime.now().hour

    print(f"Current hour: {current_hour}")
    print(f"Sunrise simulation time: {sunrise_hour}:00")
    print(f"Sunset simulation time: {sunset_hour}:00")
    print("-" * 20)

    # Check if the current time is between sunrise and sunset
    if sunrise_hour <= current_hour < sunset_hour:
        # During the day, sunlight (UV) converts the molecule to its 'cis' form.
        state = "cis"
        shape = "bent (V-shape)"
        energy = "higher-energy"
        print(f"It is daytime. The sun's UV light is present.")
        print(f"The table is in its {energy}, '{state}' form.")
        print(f"Functional State: The table is currently in a {shape}.")
    else:
        # At night, the molecule relaxes back to its more stable 'trans' form.
        state = "trans"
        shape = "straight"
        energy = "stable, lower-energy"
        print(f"It is nighttime. There is no UV light.")
        print(f"The table has relaxed to its {energy}, '{state}' form.")
        print(f"Functional State: The table has returned to a {shape} configuration.")

# Run the simulation
get_azobenzene_table_state()