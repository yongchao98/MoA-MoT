def simulate_azobenzene_table():
    """
    Simulates the daily structural change of an azobenzene-shaped picnic table.
    """

    # Properties of the two azobenzene isomers
    # The 'trans' isomer is planar and more stable. Dihedral angle is ~180 degrees.
    # The 'cis' isomer is bent and less stable. Dihedral angle is ~0 degrees.
    trans_state = {
        "name": "trans-azobenzene",
        "shape": "straight/extended",
        "dihedral_angle": 180
    }

    cis_state = {
        "name": "cis-azobenzene",
        "shape": "bent/V-shape",
        "dihedral_angle": 0
    }

    print("Simulating the daily cycle of the azobenzene picnic table:\n")

    # Sunrise: UV light from the sun causes isomerization from trans to cis
    print("Event: Sunrise (UV light exposure)")
    print("The table changes from its stable, straight form to a bent form.")
    print("Transformation Equation:")
    # Printing the change in form and the key dihedral angle number
    print(f"{trans_state['name']} (Dihedral Angle: {trans_state['dihedral_angle']}°) -> {cis_state['name']} (Dihedral Angle: {cis_state['dihedral_angle']}°)")

    print("-" * 30)

    # Sunset: In the dark, the molecule thermally relaxes back to the stable trans form
    print("Event: Sunset (Thermal relaxation in darkness)")
    print("The table slowly relaxes from its bent form back to the straight form.")
    print("Transformation Equation:")
    # Printing the change in form and the key dihedral angle number
    print(f"{cis_state['name']} (Dihedral Angle: {cis_state['dihedral_angle']}°) -> {trans_state['name']} (Dihedral Angle: {trans_state['dihedral_angle']}°)")


simulate_azobenzene_table()