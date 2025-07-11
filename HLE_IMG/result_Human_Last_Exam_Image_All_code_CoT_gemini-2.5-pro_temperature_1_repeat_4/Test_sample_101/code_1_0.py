import re

def calculate_molecular_weight():
    """
    Calculates the molecular weight of Compound A, (3-hydroxypyridin-2-yl)(phenylamino)acetonitrile.
    """
    # Standard atomic weights
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'O': 15.999
    }

    formula = "C13H11N3O"
    
    print("The final product, Compound A, is (3-hydroxypyridin-2-yl)(phenylamino)acetonitrile.")
    print(f"Its molecular formula is {formula}.")
    print("-" * 30)

    # Use regex to parse the formula and get atom counts
    # Finds pairs of (Element Symbol, Count)
    # e.g., for C13H11N3O it finds ('C', '13'), ('H', '11'), ('N', '3'), ('O', '')
    atom_counts_list = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    
    atom_counts = {}
    for element, count in atom_counts_list:
        # If count is empty, it means 1 atom
        atom_counts[element] = int(count) if count else 1

    print("The number of atoms for each element in Compound A is:")
    for element, count in atom_counts.items():
        print(f"- {element}: {count}")
    
    print("-" * 30)

    total_mw = 0
    equation_parts = []
    
    # Build the equation string and calculate the total molecular weight
    for element, count in atom_counts.items():
        weight = atomic_weights[element]
        total_mw += count * weight
        equation_parts.append(f"({count} * {weight})")
        
    final_equation = " + ".join(equation_parts)

    print("The molecular weight calculation is:")
    print(f"Molecular Weight = {final_equation}")
    print(f"\nTotal Molecular Weight of Compound A = {total_mw:.3f} g/mol")

if __name__ == '__main__':
    calculate_molecular_weight()