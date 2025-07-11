def get_product_iupac_name():
    """
    This function determines and prints the IUPAC name for the product of the
    Oxy-Cope rearrangement shown in the image.
    The final name is constructed using the determined locants (numbers).
    """
    
    # Locant for the substituent on the cyclohexanone ring
    ring_position = 3
    
    # Locants for the groups within the side chain name "4-methoxybut-3-en-2-yl"
    methoxy_position = 4
    ene_position = 3
    yl_position = 2
    
    # Construct the final IUPAC name string
    # The format is: position-(substituent)parent
    product_name = f"{ring_position}-({methoxy_position}-methoxybut-{ene_position}-en-{yl_position}-yl)cyclohexanone"
    
    print("The IUPAC name of the product is:")
    print(product_name)

get_product_iupac_name()