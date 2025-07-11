import time

def simulate_azobenzene_table():
    """
    Simulates the functional change of an azobenzene-shaped picnic table
    based on the rising and setting of the sun.
    """
    print("Simulating the daily cycle of the Azobenzene Picnic Table...\n")

    # --- SUNRISE ---
    print("Time: Sunrise")
    print("Sunlight (with UV radiation) appears.")
    print("The molecule absorbs UV light, triggering photoisomerization.")
    print("Reaction: trans-azobenzene -> cis-azobenzene")
    print("Result for the table:")
    print("The table reconfigures itself from a flat shape to a bent shape.")
    print("Current State: Bent (cis-isomer configuration)\n")

    time.sleep(1) # Pause for dramatic effect

    # --- SUNSET ---
    print("Time: Sunset")
    print("Sunlight (with UV radiation) fades.")
    print("In the dark, the molecule undergoes thermal relaxation.")
    print("Reaction: cis-azobenzene -> trans-azobenzene")
    print("Result for the table:")
    print("The table relaxes back from its bent shape to a flat shape.")
    print("Current State: Flat (trans-isomer configuration)")

# Run the simulation
simulate_azobenzene_table()