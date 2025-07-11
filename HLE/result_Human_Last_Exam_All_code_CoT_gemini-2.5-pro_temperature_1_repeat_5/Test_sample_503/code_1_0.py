import math

def solve_adjunction_count():
    """
    Calculates the number of internal adjunctions in the simplex category
    from [n] to [m].
    """
    n = 23
    m = 37

    # The number of internal adjunctions from [n] to [m] is given by the
    # binomial coefficient C(n+m, n).
    # This is equivalent to counting the number of order-preserving maps
    # f: [n] -> [m] such that f(0) = 0.
    
    total_items = n + m
    items_to_choose = n
    
    # Calculate the binomial coefficient C(total_items, items_to_choose)
    # This is C(60, 23)
    result = math.comb(total_items, items_to_choose)
    
    # Print the equation as requested
    print(f"The number of internal adjunctions from [{n}] to [{m}] is given by the binomial coefficient C({total_items}, {items_to_choose}).")
    print(f"C({total_items}, {items_to_choose}) = {total_items}! / ({items_to_choose}! * ({total_items - items_to_choose})!)")
    print(f"C({total_items}, {items_to_choose}) = {result}")

solve_adjunction_count()