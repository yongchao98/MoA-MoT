import time

class AzobenzeneTable:
    """
    A class to simulate a picnic table based on an azobenzene molecule.
    The table can change its shape based on light exposure.
    """
    def __init__(self):
        # By default, in the dark, the molecule is in its most stable, flat 'trans' form.
        self.isomer_form = "trans"
        self.shape = "flat"
        print("Before sunrise, the table is in its stable state.")
        self.print_status()

    def print_status(self):
        """Prints the current state and shape of the table."""
        print(f"  - Isomer Form: {self.isomer_form}")
        print(f"  - Table Shape: {self.shape}")

    def expose_to_sunlight(self):
        """Simulates sunrise, triggering isomerization to the 'cis' form."""
        print("\nSUNRISE: The sun's UV light shines on the table.")
        print("The molecule absorbs UV energy, causing a photoisomerization event.")
        self.isomer_form = "cis"
        self.shape = "bent"
        print("The table changes its configuration...")
        self.print_status()

    def place_in_darkness(self):
        """Simulates sunset, triggering thermal relaxation back to the 'trans' form."""
        print("\nSUNSET: The sun goes down, and the table is in darkness.")
        print("Without light energy, the molecule thermally relaxes back to its most stable form.")
        self.isomer_form = "trans"
        self.shape = "flat"
        print("The table changes its configuration...")
        self.print_status()


# --- Main Simulation ---
# Create our functional art picnic table
functional_table = AzobenzeneTable()

# Simulate a day's cycle
time.sleep(1)
functional_table.expose_to_sunlight()
time.sleep(1)
functional_table.place_in_darkness()
