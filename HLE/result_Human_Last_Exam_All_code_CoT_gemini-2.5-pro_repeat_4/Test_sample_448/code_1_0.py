import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a new string buffer
sys.stdout = captured_output = io.StringIO()

def find_rarest_noble_gas():
    """
    Identifies the rarest noble gas and calculates its percentage abundance in Earth's crust.
    The abundance data is based on mass fraction in the Earth's crust, sourced from
    "Nature's Building Blocks" by John Emsley (2001), which is suitable for the year 2002.
    """
    # Dictionary of noble gases and their mass fraction in the Earth's crust
    noble_gas_abundance_mass_fraction = {
        'Helium': 8e-9,
        'Neon': 7e-11,
        'Argon': 4e-6,
        'Krypton': 1e-10,
        'Xenon': 3e-11,
        'Radon': 4e-19
    }

    # Find the name of the gas with the minimum abundance
    rarest_gas_name = min(noble_gas_abundance_mass_fraction, key=noble_gas_abundance_mass_fraction.get)

    # Get the abundance value (mass fraction) of the rarest gas
    rarest_gas_fraction = noble_gas_abundance_mass_fraction[rarest_gas_name]

    # Convert the mass fraction to a percentage
    rarest_gas_percentage = rarest_gas_fraction * 100

    print(f"The rarest noble gas on Earth is {rarest_gas_name}.")
    print("Its approximate abundance as a percentage of terrestrial matter (crustal mass) is calculated as follows:")
    print(f"Mass Fraction: {rarest_gas_fraction}")
    print(f"Calculation: {rarest_gas_fraction} * 100")
    # Using scientific notation with one decimal place for a clean output
    print(f"Result: {rarest_gas_percentage:.1e} %")

    # The final answer is the numeric percentage value
    global final_answer
    final_answer = f"{rarest_gas_percentage:.1e}"


find_rarest_noble_gas()

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the captured output to the user
print(output)

# Print the final answer in the required format
print(f"<<<{final_answer}>>>")