def simulate_azobenzene_table():
    """
    Simulates the daily functional change of an azobenzene-shaped picnic table.
    """
    
    # Define the two states of the molecule/table
    trans_state_shape = "extended and relatively flat"
    cis_state_shape = "bent and V-shaped"
    
    print("--- Simulating the Azobenzene Picnic Table's Daily Cycle ---")
    
    # Start the simulation at night, in the most stable state
    print("\nInitial State (During the night):")
    print(f"The table is in its stable trans-isomer form: {trans_state_shape}.")

    # --- Sunrise Event ---
    print("\nSUNRISE:")
    print("The sun rises, emitting UV light.")
    print("The azobenzene molecule absorbs the UV photons.")
    print("This triggers a change from the trans to the cis isomer.")
    
    # Print the "equation" for this photochemical reaction
    # The final equation is: trans-azobenzene --(UV light)--> cis-azobenzene
    print("\nThe functional change is described by the equation:")
    print("trans-azobenzene --> cis-azobenzene")
    print("Trigger: UV light")
    
    print(f"\nResult: The table physically reconfigures itself from a '{trans_state_shape}' shape to a '{cis_state_shape}' shape.")

    # --- Sunset Event ---
    print("\nSUNSET:")
    print("The sun sets, and the source of UV light disappears.")
    print("The less stable cis-isomer relaxes back to the stable trans-isomer.")
    print("This change is driven by thermal energy (heat) and ambient visible light.")
    
    # Print the "equation" for this reverse reaction
    # The final equation is: cis-azobenzene --(Heat / Visible Light)--> trans-azobenzene
    print("\nThe functional change is described by the equation:")
    print("cis-azobenzene --> trans-azobenzene")
    print("Trigger: Heat / Visible Light")

    print(f"\nResult: The table physically reconfigures itself from a '{cis_state_shape}' shape back to a '{trans_state_shape}' shape for the night.")


simulate_azobenzene_table()