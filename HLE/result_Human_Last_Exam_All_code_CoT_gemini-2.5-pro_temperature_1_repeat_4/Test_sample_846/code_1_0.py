def simulate_azobenzene_table(time_of_day):
    """
    Simulates the state of an azobenzene-shaped table based on the time of day.
    """
    # Trans-azobenzene is planar (flat) and more stable.
    # Cis-azobenzene is bent and less stable.
    trans_state = "trans (flat)"
    cis_state = "cis (bent)"

    print(f"--- It is now {time_of_day}. ---")

    if time_of_day.lower() == "sunrise":
        print("The sun is rising, exposing the table to UV light.")
        print(f"The molecule-table is undergoing photoisomerization.")
        print(f"State change: {trans_state} -> {cis_state}")
        print("The table should physically change from a flat to a bent shape.")

    elif time_of_day.lower() == "daytime":
        print("The sun is out, and UV light is present.")
        print(f"The table remains in the {cis_state} state.")

    elif time_of_day.lower() == "sunset":
        print("The sun is setting, and UV light is fading.")
        print("The molecule-table is undergoing thermal relaxation.")
        print(f"State change: {cis_state} -> {trans_state}")
        print("The table should physically relax from a bent back to a flat shape for the night.")
        
    elif time_of_day.lower() == "night":
        print("It is dark.")
        print(f"The table remains in its stable {trans_state} state.")

    else:
        print("Unknown time of day.")

# Simulate a full day cycle
simulate_azobenzene_table("sunrise")
print("\n")
simulate_azobenzene_table("daytime")
print("\n")
simulate_azobenzene_table("sunset")
print("\n")
simulate_azobenzene_table("night")