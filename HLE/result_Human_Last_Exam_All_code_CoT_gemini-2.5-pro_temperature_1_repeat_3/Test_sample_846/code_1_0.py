def simulate_azobenzene_table(time_of_day):
    """
    Simulates the behavior of an azobenzene-shaped table
    based on the time of day and corresponding light conditions.
    """
    if time_of_day.lower() == "sunrise":
        print("Event: The sun rises.")
        print("Trigger: The table is exposed to UV radiation from sunlight.")
        print("Action: To be functionally accurate, the molecule must undergo photoisomerization.")
        print("Result: The table changes shape from its stable, elongated 'trans' form to its bent 'cis' form.")
        print("-" * 30)
    elif time_of_day.lower() == "sunset":
        print("Event: The sun sets.")
        print("Trigger: The UV radiation source is removed.")
        print("Action: The molecule must thermally relax back to its most stable state.")
        print("Result: The table changes shape from its bent 'cis' form back to the elongated 'trans' form.")
        print("-" * 30)
    else:
        print(f"At {time_of_day}, the table remains in its current state until the next sunrise or sunset.")

# Demonstrate the daily cycle
simulate_azobenzene_table("sunrise")
simulate_azobenzene_table("sunset")
