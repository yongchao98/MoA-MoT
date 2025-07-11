def solve_century_calculation():
    """
    Calculates the century for a given year and prints the process.
    The medallion is dated to the reign of Roman Emperor Galerius (305-311 AD).
    We will use a representative year, 310 AD, for the calculation.
    """
    year = 310
    
    # The formula to calculate the century from a year is (year - 1) // 100 + 1
    # Step 1: Subtract 1 from the year
    step1_val = year - 1
    
    # Step 2: Integer divide by 100
    step2_val = step1_val // 100
    
    # Step 3: Add 1
    century = step2_val + 1

    print(f"The medallion dates to approximately {year} A.D.")
    print(f"To find the century, we calculate: ({year} - 1) // 100 + 1")
    print(f"Which is: ({step1_val}) // 100 + 1")
    print(f"Which is: {step2_val} + 1")
    print(f"Resulting in the {century}th century A.D.")

solve_century_calculation()