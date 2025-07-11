import time

class AzobenzeneTable:
    """
    A class to simulate the behavior of a functional azobenzene picnic table.
    """
    def __init__(self):
        # The table starts in its most stable, default state (trans).
        self.state = "trans (extended/straight)"
        print(f"The table is constructed. Initial state: {self.state}")

    def sunrise_event(self):
        """Simulates the effect of sunrise (UV light exposure)."""
        print("\nThe sun rises, exposing the table to UV light...")
        print("This triggers photoisomerization.")
        old_state = self.state
        self.state = "cis (bent/folded)"
        print(f"The table's shape changes from {old_state} to {self.state}.")

    def sunset_event(self):
        """Simulates the effect of sunset (darkness and thermal relaxation)."""
        print("\nThe sun sets, and darkness falls...")
        print("Thermal relaxation begins.")
        old_state = self.state
        self.state = "trans (extended/straight)"
        print(f"The table's shape relaxes from {old_state} back to {self.state}.")

# --- Main execution to demonstrate the daily cycle ---
if __name__ == "__main__":
    functional_table = AzobenzeneTable()
    
    # Simulate the day's events
    functional_table.sunrise_event()
    functional_table.sunset_event()
