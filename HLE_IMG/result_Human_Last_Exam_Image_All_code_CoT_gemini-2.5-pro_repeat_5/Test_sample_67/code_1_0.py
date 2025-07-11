def solve_electron_energy():
    """
    Calculates and prints the minimum energy of electron 1 for the described process.

    The derivation from conservation laws under threshold conditions leads to:
    E1_min = Eg + (3/2)*Eg
    """

    # We represent the band gap energy symbolically. The coefficients are numbers.
    coefficient_rest_energy = 1
    coefficient_kinetic_energy = 3/2

    # The total energy is the sum of the rest energy (at k=0) and the kinetic energy.
    total_coefficient = coefficient_rest_energy + coefficient_kinetic_energy

    print("The minimum energy for electron 1 (E1_min) is the sum of its energy at the bottom of the band (Eg) and the minimum required kinetic energy.")
    print("From the laws of conservation of energy and momentum, the minimum kinetic energy is found to be (3/2)*Eg.")
    print("\nSo, the final equation is:")
    print(f"E1_min = {coefficient_rest_energy}*Eg + ({coefficient_kinetic_energy})*Eg")
    print(f"E1_min = ({coefficient_rest_energy + coefficient_kinetic_energy})*Eg")
    print(f"\nTherefore, the minimum energy of electron 1 must be {total_coefficient} times the band gap energy, Eg.")

solve_electron_energy()