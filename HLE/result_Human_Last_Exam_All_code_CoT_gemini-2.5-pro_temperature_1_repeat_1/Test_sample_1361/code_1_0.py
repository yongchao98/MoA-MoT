def solve_ballet_riddle():
    """
    This function calculates the number of gargouillades in Darcey Bussell's
    2005 performance of Sylvia's Act III solo.

    The calculation is based on the years provided in the prompt.
    """
    year_created = 1952
    year_performed = 2005

    # Step 1: Calculate the difference between the performance year and the creation year.
    year_difference = year_performed - year_created

    # Step 2: Sum the digits of the year difference.
    # Convert the number to a string to access its digits.
    s_year_difference = str(year_difference)
    
    # Initialize a list to hold the digits for the equation display
    digits = [int(d) for d in s_year_difference]
    
    # Calculate the total
    num_gargouillades = sum(digits)

    # Step 3: Print the equation as requested.
    equation_parts = [str(d) for d in digits]
    equation_str = " + ".join(equation_parts)
    
    print(f"The number of gargouillades is derived from the years:")
    print(f"Years between performance and creation: {year_performed} - {year_created} = {year_difference}")
    print(f"Summing the digits of the difference: {equation_str} = {num_gargouillades}")

solve_ballet_riddle()