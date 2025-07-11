def solve_silicate_complexity():
    """
    Identifies mineral specimens, classifies their silicate structure,
    and orders them by increasing structural complexity.
    """
    # Step 1 & 2: Identify specimens and classify their silicate structure.
    # A complexity_rank is assigned for sorting purposes (lower is less complex).
    specimens = [
        {'id': 'A', 'mineral': 'Pyroxene', 'structure': 'Inosilicate (Single Chain)', 'complexity_rank': 2},
        {'id': 'B', 'mineral': 'Mica', 'structure': 'Phyllosilicate (Sheet)', 'complexity_rank': 3},
        {'id': 'C', 'mineral': 'Olivine', 'structure': 'Nesosilicate (Isolated Tetrahedra)', 'complexity_rank': 1},
        {'id': 'D', 'mineral': 'Quartz', 'structure': 'Tectosilicate (Framework)', 'complexity_rank': 4}
    ]

    print("Step 1: Identifying each specimen and its silicate structure.")
    for specimen in specimens:
        print(f"Specimen {specimen['id']}: Identified as {specimen['mineral']}, a {specimen['structure']}.")
    
    print("\nStep 2: Ordering the specimens by increasing silicate structure complexity.")
    print("Complexity Order: Nesosilicate < Inosilicate < Phyllosilicate < Tectosilicate")

    # Step 3: Sort the specimens based on their complexity rank.
    sorted_specimens = sorted(specimens, key=lambda x: x['complexity_rank'])
    
    # Extract the IDs for the final answer
    ordered_ids = [s['id'] for s in sorted_specimens]
    
    print("\nFinal Order:")
    print("The specimens ordered by increasing silicate structure complexity are:")
    # Using ' -> '.join() to create the arrow-separated string
    print(' -> '.join(ordered_ids))

# Run the function to display the solution
solve_silicate_complexity()