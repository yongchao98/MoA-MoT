import sys
# This script explains why temperature is the essential parameter.
# In statistical mechanics, the average kinetic energy of a gas molecule is
# directly proportional to the absolute temperature (T).
# The formula is: Average Kinetic Energy = (3/2) * k * T
# where 'k' is the Boltzmann constant.

def explain_role_of_temperature(temperature_kelvin):
    """
    Explains the state of gas molecules at a given temperature.
    """
    if temperature_kelvin < 0:
        print("Error: Temperature cannot be below absolute zero.", file=sys.stderr)
        return

    print(f"Analyzing the system at a temperature of {temperature_kelvin} Kelvin...")

    if temperature_kelvin == 0:
        print("At absolute zero, the average kinetic energy of gas molecules is zero.")
        print("This means the molecules are effectively motionless.")
        print("Without motion, they cannot pass through the one-way door.")
        print("Therefore, the gas cannot be trapped on one side.")
    else:
        print(f"At a non-zero temperature of {temperature_kelvin} K, molecules have kinetic energy.")
        print("This energy causes them to be in constant, random motion.")
        print("This random motion will eventually lead to all molecules passing through the one-way door,")
        print("trapping the gas entirely in one compartment.")

# --- Main Demonstration ---
print("Demonstrating the effect of temperature on molecular motion:")
print("="*60)

# Case 1: Absolute Zero
explain_role_of_temperature(0)

print("\n" + "="*60 + "\n")

# Case 2: A non-zero temperature (e.g., room temperature)
explain_role_of_temperature(298)

print("\n" + "="*60)
print("\nConclusion: A non-zero Temperature is the required experimental parameter to provide the kinetic energy for molecular motion.")
