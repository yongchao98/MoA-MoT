import sys

def propose_methanol_charges():
    """
    Proposes and verifies a reasonable set of partial charges for methanol.
    """
    # Choice D: A chemically intuitive, symmetric, and charge-neutral set.
    charges = {
        "Carbon": 0.1450,
        "Methyl H1": 0.0400,
        "Methyl H2": 0.0400,
        "Methyl H3": 0.0400,
        "Oxygen": -0.6830,
        "Hydroxyl hydrogen": 0.4180
    }

    print("Proposed partial charge assignments for a new methanol model:")
    print("-" * 60)
    # Use a consistent field width for better alignment
    for atom, charge in charges.items():
        print(f"{atom:<20}: {charge:>8.4f}")
    print("-" * 60)

    # Verify charge neutrality and print the equation
    total_charge = sum(charges.values())
    
    # Using sys.stdout.write to avoid spaces between print arguments
    print("\nVerification of charge neutrality:")
    charge_values = list(charges.values())
    
    # Print the equation with each number explicitly
    sys.stdout.write(f"{charge_values[0]:.4f} (C) + ")
    sys.stdout.write(f"{charge_values[1]:.4f} (H1) + ")
    sys.stdout.write(f"{charge_values[2]:.4f} (H2) + ")
    sys.stdout.write(f"{charge_values[3]:.4f} (H3) + ")
    sys.stdout.write(f"{charge_values[4]:.4f} (O) + ")
    sys.stdout.write(f"{charge_values[5]:.4f} (H_ox) = ")
    print(f"{total_charge:.4f}")

    # Check if the total charge is effectively zero
    if abs(total_charge) < 1e-9:
        print("\nThe total charge is zero. This is a valid and reasonable set of charges.")
    else:
        print(f"\nWarning: The total charge is {total_charge:.4f}, not zero.")

if __name__ == '__main__':
    propose_methanol_charges()
