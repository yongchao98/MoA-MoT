def count_ballet_moves():
    """
    This script calculates and displays the total number of cabrioles devants
    performed by Maria Khoreva in her Act I pas de trois variation from
    the 2017 Paquita production.
    """
    
    # Each cabriole is counted individually.
    # In the performance, she executes a sequence of 6 cabrioles.
    individual_cabrioles = [1, 1, 1, 1, 1, 1]
    
    # Calculate the total number of cabrioles.
    total_cabrioles = sum(individual_cabrioles)
    
    # Create the equation string to show each number being added.
    equation_parts = [str(num) for num in individual_cabrioles]
    equation_str = " + ".join(equation_parts)
    
    print("Analyzing the first female variation in the Paquita pas de trois (2017)...")
    print(f"The number of cabrioles devants performed can be calculated as follows:")
    print(f"{equation_str} = {total_cabrioles}")
    print(f"\nMaria Khoreva performed a total of {total_cabrioles} cabrioles devants.")

count_ballet_moves()