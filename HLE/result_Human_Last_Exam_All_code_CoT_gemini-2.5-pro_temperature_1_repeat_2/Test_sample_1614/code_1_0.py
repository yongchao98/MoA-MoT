def solve_tiling():
    """
    Calculates the number of ways to tile a 2x4 board, T(4).
    """
    # The recurrence relation is T(n) = T(n-1) + T(n-2) + T(n-4).
    # Base cases: T(0) = 1, T(n) = 0 for n < 0.
    
    # T(0)
    t0 = 1
    
    # T(1) = T(0) + T(-1) + T(-3) = 1 + 0 + 0 = 1
    t1 = t0 
    
    # T(2) = T(1) + T(0) + T(-2) = 1 + 1 + 0 = 2
    t2 = t1 + t0
    
    # T(3) = T(2) + T(1) + T(-1) = 2 + 1 + 0 = 3
    t3 = t2 + t1
    
    # T(4) = T(3) + T(2) + T(0) = 3 + 2 + 1 = 6
    t4 = t3 + t2 + t0
    
    print(f"The number of ways to tile a 2x4 board is T(4).")
    print(f"The recurrence relation is T(n) = T(n-1) + T(n-2) + T(n-4).")
    print(f"We calculate step-by-step:")
    print(f"T(0) = 1")
    print(f"T(1) = T(0) = 1")
    print(f"T(2) = T(1) + T(0) = {t1} + {t0} = {t2}")
    print(f"T(3) = T(2) + T(1) = {t2} + {t1} = {t3}")
    print(f"T(4) = T(3) + T(2) + T(0) = {t3} + {t2} + {t0} = {t4}")

solve_tiling()