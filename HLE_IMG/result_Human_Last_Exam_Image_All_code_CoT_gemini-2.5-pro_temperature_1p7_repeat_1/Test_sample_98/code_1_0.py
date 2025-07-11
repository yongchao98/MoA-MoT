def find_century():
    """
    This function determines and prints the century of the Kerch Medallion.
    The medallion is dated to the reign of Constantius II, specifically around 350-360 A.D.
    Years from 301 to 400 A.D. fall within the 4th century A.D.
    """
    century = 4
    era = "A.D."
    
    print(f"The medallion is ascribed to the {century}th century {era}.")
    # To match the requested output format like "9 A.D."
    print(f"{century} {era}")

find_century()