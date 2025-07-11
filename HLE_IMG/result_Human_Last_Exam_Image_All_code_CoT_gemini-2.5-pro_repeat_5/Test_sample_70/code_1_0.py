def identify_pericyclic_reactions():
    """
    Identifies and describes the two pericyclic reactions in the given thermal transformation.
    The transformation proceeds from a bicyclo[6.2.0]decapentaene to 9,10-dihydronaphthalene.
    """
    
    # The first reaction is the opening of the four-membered ring.
    # This involves 1 pi bond and 1 sigma bond, for a total of 4 electrons.
    reaction1_electrons = 4
    reaction1_type = "electrocyclic ring opening"
    
    # The second reaction is the closure of the ten-membered ring intermediate.
    # This involves a hexatriene system with 3 pi bonds, for a total of 6 electrons.
    reaction2_electrons = 6
    reaction2_type = "electrocyclic ring closure"
    
    print("The thermal transformation involves two sequential pericyclic reactions:")
    print(f"1. A [{reaction1_electrons}π] {reaction1_type}.")
    print(f"2. A [{reaction2_electrons}π] {reaction2_type}.")

if __name__ == '__main__':
    identify_pericyclic_reactions()