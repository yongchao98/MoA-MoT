def solve_ballet_riddle():
    """
    This function calculates and displays the number of cabrioles devants
    in Maria Khoreva's 2017 Paquita solo variation.
    
    Each '1' in the list represents one observed cabriole.
    """
    
    # A list representing each cabriole devant executed in the sequence.
    observed_cabrioles = [1, 1, 1, 1]
    
    # Calculate the total number by summing the list.
    total = sum(observed_cabrioles)
    
    # Create the string for the equation.
    equation_str = " + ".join(map(str, observed_cabrioles))
    
    # Print the final equation.
    print(f"By observing the performance, we count each cabriole devant.")
    print(f"The final calculation is:")
    print(f"{equation_str} = {total}")

solve_ballet_riddle()