import networkx as nx

def check_cope_rearrangement_product():
    """
    Checks the correctness of the product from the aza-Cope rearrangement of
    (1S,4R)-2-vinyl-2-azabicyclo[2.2.1]hept-5-ene.
    """

    # Step 1: Define the reactant molecule as a graph.
    # Atom labels are based on IUPAC numbering for the bicyclo[2.2.1]heptane core,
    # with C8 and C9 representing the vinyl group (-C8H=C9H2) on N2.
    reactant = nx.Graph()
    reactant_bonds = [
        ('C1', 'N2'), ('C1', 'C6'), ('C1', 'C7'), ('C1', 'C4'),
        ('N2', 'C3'), ('N2', 'C8'),
        ('C3', 'C4'),
        ('C4', 'C5'), ('C4', 'C7'),
        ('C5', 'C6', {'order': 2}),
        ('C8', 'C9', {'order': 2})
    ]
    reactant.add_edges_from(reactant_bonds)

    # Step 2: Identify the 1,5-diene system and apply the [3,3] shift.
    # The system is C9=C8-N2-C1-C6=C5.
    # Let's create the predicted product by transforming the reactant graph.
    predicted_product = reactant.copy()

    # Rule 1: Break the bond between atoms 3 and 4 (N2-C1).
    predicted_product.remove_edge('N2', 'C1')

    # Rule 2: Form a bond between atoms 1 and 6 (C9-C5).
    predicted_product.add_edge('C9', 'C5', order=1)

    # Rule 3: Shift the double bonds.
    # C8=C9 becomes C8-C9, and N2-C8 becomes N2=C8.
    predicted_product['C8']['C9']['order'] = 1
    predicted_product.add_edge('C8', 'N2', order=2)
    # C5=C6 becomes C5-C6, and C1-C6 becomes C1=C6.
    predicted_product['C5']['C6']['order'] = 1
    predicted_product.add_edge('C1', 'C6', order=2)

    # Step 3: Analyze the predicted product's skeleton.
    # The product is a bicyclo[4.3.0]nonane derivative (a fused 5- and 6-membered ring).
    # The fusion carbons are C4 and C5.
    # The nitrogen atom (N2) is in the 6-membered ring.
    # By IUPAC rules for bicyclic systems, this structure is a 3-azabicyclo[4.3.0]nona-3,8-diene.

    # Step 4: Define the structure from the given answer D.
    # Answer D: 4,4a,5,7a-tetrahydro-3H-cyclopenta[c]pyridine.
    # This name corresponds to a 2-azabicyclo[4.3.0]nonane skeleton.
    # The double bonds are at positions C1=N2 and C6=C7 (using the fused system's IUPAC numbering).
    answer_d = nx.Graph()
    answer_d_bonds = [
        ('C1_ans', 'N2_ans', {'order': 2}), ('C1_ans', 'C7a_ans'),
        ('N2_ans', 'C3_ans'),
        ('C3_ans', 'C4_ans'),
        ('C4_ans', 'C4a_ans'),
        ('C4a_ans', 'C5_ans'), ('C4a_ans', 'C7a_ans'),
        ('C5_ans', 'C6_ans'),
        ('C6_ans', 'C7_ans', {'order': 2}),
        ('C7_ans', 'C7a_ans')
    ]
    answer_d.add_edges_from(answer_d_bonds)

    # Step 5: Compare the predicted structure with the answer's structure.
    # The predicted product is a 3-aza isomer, while the answer is a 2-aza isomer.
    # We can confirm they are different by checking the neighborhood of the nitrogen atom.
    # In the predicted product, the nitrogen (N2) is bonded to C3 and C8.
    # In the answer's structure, the nitrogen (N2_ans) is bonded to C1_ans and C3_ans.
    # These describe different connectivities. They are constitutional isomers, not the same molecule.

    # A graph isomorphism check would formally prove they are different, but the difference
    # in the nitrogen's position in the bicyclic skeleton is sufficient.
    
    is_isomorphic = nx.is_isomorphic(predicted_product, answer_d, node_match=lambda n1, n2: True)

    if is_isomorphic:
        return "Correct"
    else:
        return (
            "The provided answer is incorrect based on the standard mechanism for a Cope rearrangement.\n"
            "1. The aza-Cope rearrangement is a [3,3]-sigmatropic reaction. Applying this mechanism to the reactant correctly yields a product with a `bicyclo[4.3.0]nonane` skeleton (a fused 5- and 6-membered ring system).\n"
            "2. However, the specific connectivity (i.e., the position of the nitrogen and the double bonds) predicted by the mechanism does not match the structure described by the name in answer D.\n"
            "3. The mechanism predicts a `3-azabicyclo[4.3.0]nona-3,8-diene` isomer.\n"
            "4. Answer D, `4,4a,5,7a-tetrahydro-3H-cyclopenta[c]pyridine`, corresponds to a `2-azabicyclo[4.3.0]nona-1,6-diene` isomer.\n"
            "Since the predicted product and the given answer are different constitutional isomers, the answer does not satisfy the constraints of the specified reaction mechanism."
        )

result = check_cope_rearrangement_product()
print(result)