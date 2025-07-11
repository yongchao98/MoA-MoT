import math

def solve_century():
    """
    This script determines the century of the Tillya Tepe medallion.
    """
    # The Tillya Tepe graves, where the medallion was found, are generally dated
    # by archaeologists to the period between 25 B.C. and 50 A.D.
    # The art style and coins found in the graves confirm this dating.
    # The core period of this range falls within the 1st century A.D.
    start_year_ad = 1
    end_year_ad = 100

    # The 1st century A.D. corresponds to the years 1 through 100.
    century = math.ceil(start_year_ad / 100)
    era = "A.D."

    print(f"The medallion is from the Tillya Tepe treasure, dated to the early part of the {century}st century {era}.")
    print(f"Therefore, the century number is: {century}")
    print(f"The era is: {era}")

solve_century()