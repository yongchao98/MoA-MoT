def generate_lc_molecule(n):
    """
    Generates the name and formula for a potential liquid crystal molecule
    based on the design CnH2n+1-Ph-CN.
    
    Args:
        n (int): The number of carbon atoms in the alkyl chain.
    
    Returns:
        tuple: A tuple containing the chemical name and the formula string.
    """
    # A dictionary to map chain length to the standard chemical prefix
    prefixes = {
        1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl", 
        5: "pentyl", 6: "hexyl", 7: "heptyl", 8: "octyl"
    }
    
    alkyl_name = prefixes.get(n, f"C{n}-alkyl")
    # The structure is a para-substituted (1,4) benzonitrile
    chemical_name = f"4-{alkyl_name}benzonitrile"
    
    # Calculate the number of hydrogen atoms in the alkyl chain (CnH2n+1)
    hydrogen_count = 2 * n + 1
    
    # The final equation representing the molecule's structure
    # Where Ph represents the C6H4 benzene ring
    formula = f"C({n})H({hydrogen_count})-Ph-CN"
    
    return chemical_name, formula

print("--- Liquid Crystal Design ---")
print("Goal: Single benzene ring molecule with transition temperature near room temperature.")
print("General Structure: C(n)H(2n+1)-Ph-CN")
print("Tuning Method: Vary the alkyl chain length 'n'.")
print("-" * 30)

# As per the design guide, n=5 is a good starting point for room temperature behavior.
n_start = 5
name, formula = generate_lc_molecule(n_start)
print(f"Suggested Starting Point (n = {n_start}):")
print(f"  Name: {name}")
print(f"  Final Equation: {formula}")
print("  This molecule is a prime candidate for room temperature liquid crystal phases.")
print("-" * 30)

# To raise the transition temperature, shorten the chain.
n_short = 4
name, formula = generate_lc_molecule(n_short)
print(f"To Increase Transition Temp (n = {n_short}):")
print(f"  Name: {name}")
print(f"  Final Equation: {formula}")
print("-" * 30)

# To lower the transition temperature, lengthen the chain.
n_long = 6
name, formula = generate_lc_molecule(n_long)
print(f"To Decrease Transition Temp (n = {n_long}):")
print(f"  Name: {name}")
print(f"  Final Equation: {formula}")
print("-" * 30)
