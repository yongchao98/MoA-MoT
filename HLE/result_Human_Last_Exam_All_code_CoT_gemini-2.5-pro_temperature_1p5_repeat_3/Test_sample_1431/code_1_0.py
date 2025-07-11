def calculate_monk_force():
    """
    Calculates the force F required to lift a rope of mass m and length l
    such that it has a speed v when it just leaves the ground.
    
    The user should modify the placeholder values for mass, length, and speed below.
    """
    # --- Please modify these values for your specific rope ---
    m = 10.0  # mass of the rope in kilograms (kg)
    l = 20.0  # length of the rope in meters (m)
    v = 5.0   # final speed in meters per second (m/s)
    # ---------------------------------------------------------

    g = 9.8   # acceleration due to gravity in m/s^2

    # Calculate the components of the Work-Energy equation
    delta_potential_energy = m * g * (l / 2)
    delta_kinetic_energy = 0.5 * m * v**2
    total_work = delta_potential_energy + delta_kinetic_energy
    
    # Calculate the final force F = Work / distance
    force = total_work / l

    # Print the explanation and the final equation with values
    print("The solution is found using the Work-Energy Theorem: W_force = ΔPotential_Energy + ΔKinetic_Energy")
    print("\n1. Work done by the applied force F over a distance l:")
    print(f"   W_force = F * l")
    
    print("\n2. Change in Potential Energy (ΔU):")
    print(f"   ΔU = m * g * (l/2) = {m} * {g} * ({l}/2) = {delta_potential_energy:.2f} Joules")

    print("\n3. Change in Kinetic Energy (ΔK):")
    print(f"   ΔK = (1/2) * m * v^2 = 0.5 * {m} * {v}^2 = {delta_kinetic_energy:.2f} Joules")

    print("\nPutting it all together:")
    print(f"   F * l = ΔU + ΔK")
    print(f"   F * {l} = {delta_potential_energy:.2f} + {delta_kinetic_energy:.2f}")
    print(f"   F * {l} = {total_work:.2f}")
    print(f"   F = {total_work:.2f} / {l}")

    print("\n---------------------------------------------------")
    print(f"The exact mystical force F required is: {force:.2f} Newtons")
    print("---------------------------------------------------")

# Execute the calculation
calculate_monk_force()