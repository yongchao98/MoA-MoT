def solve_sequence():
    """
    This function determines and prints the first 4 terms of the sequence D_g.
    
    D_g is the smallest degree of a finite etale cover of the moduli stack A_g of
    principally polarised abelian varieties (PPAVs) of dimension g, such that the
    universal polarisation becomes symmetric.
    
    - For g=1, the universal polarisation on A_1 is symmetric, so D_1 = 1.
    - For g>=2, the universal polarisation is not symmetric, and the obstruction is a
      non-trivial line bundle of order 2. Thus D_g = 2 for g>=2.
    """
    
    # The first four terms of the sequence D_g for g = 1, 2, 3, 4.
    D1 = 1
    D2 = 2
    D3 = 2
    D4 = 2
    
    print("The first four terms of the sequence D_g are:")
    print(f"D_1 = {D1}")
    print(f"D_2 = {D2}")
    print(f"D_3 = {D3}")
    print(f"D_4 = {D4}")

solve_sequence()