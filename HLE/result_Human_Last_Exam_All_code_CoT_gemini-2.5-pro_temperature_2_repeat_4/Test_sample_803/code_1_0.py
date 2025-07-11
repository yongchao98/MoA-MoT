import subprocess
import re

def is_prime(n):
    """A simple primality test function."""
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_groups_with_gap(order):
    """
    Calls the GAP software as a subprocess to get a list of all
    nonabelian groups of a given order from its SmallGroups library.
    """
    # This GAP command finds all nonabelian groups of the given order and
    # prints their GAP ID and structure description.
    command = (
        f"gap -q -c 'L:=AllSmallGroups({order}, IsAbelian, false);; "
        f"for g in L do Print(IdGroup(g), \";\", StructureDescription(Group(g)), \"\\n\"); od; QUIT;'"
    )
    
    try:
        # We execute the command in a shell.
        # This requires 'gap' to be in the system's PATH.
        result = subprocess.run(
            command, 
            shell=True, 
            capture_output=True, 
            text=True, 
            timeout=60,  # A timeout to prevent very long computations
            check=True   # Raise an exception if GAP returns a non-zero exit code
        )
        # Return the output, split into lines for each group
        output = result.stdout.strip()
        if not output:
            return []
        return output.split('\n')

    except FileNotFoundError:
        return ["Error: The 'gap' command was not found. Please ensure GAP is installed and that its 'bin' directory is in your system's PATH environment variable."]
    except subprocess.CalledProcessError as e:
        # This often happens if the order is too large or not in the SmallGroups library.
        return [f"Error: GAP encountered a problem. The order {order} may be too large or not available in GAP's SmallGroups library. Details: {e.stderr.strip()}"]
    except subprocess.TimeoutExpired:
        return [f"Error: The GAP process timed out. The calculation for order {order} is too complex or is taking too long."]

def main():
    """
    Main function to find nonabelian filled groups of order 2*q^m.
    """
    print("This script finds all nonabelian filled groups of order 2*q^m.")
    print("Based on established theorems, all solvable groups are filled, and groups of order 2*q^m are always solvable.")
    print("Therefore, this script simply finds all nonabelian groups of that order.")
    print("-" * 70)
    
    try:
        q_str = input("Enter an odd prime number q (e.g., 3, 5, 7): ")
        q = int(q_str)
        if q % 2 == 0 or not is_prime(q):
            print(f"Error: {q} is not an odd prime number.")
            return

        m_str = input("Enter a natural number m (>= 1): ")
        m = int(m_str)
        if m < 1:
            print("Error: m must be a natural number (1, 2, 3, ...).")
            return
            
    except ValueError:
        print("Invalid input. Please enter integers for q and m.")
        return

    # Calculate the order and display the equation
    order = 2 * (q**m)
    print(f"\nThe equation for the group order is: 2 * {q}^{m} = {order}")
    print(f"Searching for nonabelian groups of order {order} using GAP...")
    print("(This may take a moment, especially for larger orders...)\n")
    
    results = find_groups_with_gap(order)
    
    if not results:
        print(f"No nonabelian groups of order {order} were found.")
        print("This could be because no such groups exist, or the order is too large for the SmallGroups library.")
    elif "Error:" in results[0]:
        # Print the detailed error message from the helper function
        print(results[0])
    else:
        print(f"Found {len(results)} nonabelian filled group(s) of order {order}:")
        for line in results:
            # Expected line format: "[order,id];StructureDescription"
            parts = line.split(';')
            if len(parts) == 2:
                group_id_str, struct_desc = parts
                print(f"  - GAP ID {group_id_str.strip()}: {struct_desc.strip()}")
            else:
                print(f"  - {line.strip()}") # Fallback for unexpected format

if __name__ == '__main__':
    main()
