import math

def solve_medallion_date():
    """
    This function solves for the date of the medallion based on historical context.
    The medallion's findspot in Pannonia and the hint '9 A.D.' point to the
    end of the Great Illyrian Revolt.
    """
    
    # The Great Illyrian/Pannonian Revolt started in 6 A.D.
    start_year_ad = 6
    
    # The revolt lasted for approximately 3 full years.
    duration_of_revolt = 3
    
    # We calculate the year the revolt ended, which is the likely date of the medallion.
    end_year_ad = start_year_ad + duration_of_revolt
    
    print("The medallion likely commemorates the end of the Great Illyrian Revolt in Pannonia.")
    print("The revolt started in the year 6 A.D. and lasted for 3 years.")
    print("The equation for the final year is:")
    print(f"{start_year_ad} + {duration_of_revolt} = {end_year_ad}")
    print("\nTherefore, archaeologists would likely ascribe it to the year:")
    print(f"{end_year_ad} A.D.")

solve_medallion_date()