def simulate_azobenzene_table_function():
    """
    This function describes the behavior of a functional azobenzene
    picnic table throughout a day-night cycle.
    """
    # At sunrise, the table is in its stable 'trans' state.
    print("Event: Sunrise")
    print("Sunlight (containing UV rays) hits the table.")
    
    # The isomerization process is governed by light wavelength.
    # The 'equation' for this change involves absorbing a photon of a specific energy.
    trans_to_cis_wavelength = 365 # Wavelength in nanometers
    print(f"The molecule absorbs UV light, with a characteristic wavelength of ~{trans_to_cis_wavelength} nm.")
    print("This triggers 'trans-to-cis' isomerization.")
    print("Result: The flat picnic table BENDS and CONTRACTS into its 'cis' form.\n")

    # At sunset, the UV light source is removed.
    print("Event: Sunset")
    print("It is now dark, and the table is in its 'cis' form.")

    # The relaxation back to the 'trans' form can be driven by heat (darkness) or visible light.
    cis_to_trans_wavelength = 440 # Wavelength in nanometers for light-driven relaxation
    print(f"In the dark, the molecule thermally relaxes. This can also be driven by visible light of ~{cis_to_trans_wavelength} nm.")
    print("This triggers 'cis-to-trans' isomerization.")
    print("Result: The bent picnic table RELAXES and FLATTENS back to its stable 'trans' form.")

# Run the simulation
simulate_azobenzene_table_function()