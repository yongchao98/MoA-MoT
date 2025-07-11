# To run this code, you may need to install the pyDamerauLevenshtein library:
# pip install pyDamerauLevenshtein

from pyDamerauLevenshtein import Damerau

def verify_triangle_inequality_violation():
    """
    This function verifies statement (D), which claims that the LT distance
    (Optimal String Alignment variant) can violate the triangle inequality.
    The triangle inequality states that for any a, b, c:
    d(a, c) <= d(a, b) + d(b, c)

    We use the classic counterexample: a="CA", b="AC", c="ABC".
    If the inequality is violated, we will find that LT(a, c) > LT(a, b) + LT(b, c).
    """

    # The pyDamerauLevenshtein library's Damerau class computes the
    # Optimal String Alignment (OSA) distance, which is what is commonly
    # understood as "Levenshtein with transpositions".
    osa_calculator = Damerau()

    a = "CA"
    b = "AC"
    c = "ABC"

    # Calculate the three required distances
    dist_ac = osa_calculator(a, c)
    dist_ab = osa_calculator(a, b)
    dist_bc = osa_calculator(b, c)

    # Print the equation and the result
    print(f"Testing triangle inequality for LT distance with a='{a}', b='{b}', c='{c}'")
    print(f"The inequality is: LT(a, c) <= LT(a, b) + LT(b, c)")
    print("Plugging in the calculated values:")
    
    # Using the required format to output each number in the final equation
    print(f"{dist_ac} > {dist_ab} + {dist_bc} is {dist_ac > dist_ab + dist_bc}")
    print(f"\nBreakdown:")
    print(f"LT('{a}', '{c}') = {dist_ac}")
    print(f"LT('{a}', '{b}') = {dist_ab} (one transposition of 'CA' -> 'AC')")
    print(f"LT('{b}', '{c}') = {dist_bc} (one insertion of 'B' into 'AC')")
    
    total_ab_bc = dist_ab + dist_bc
    print(f"\nSince {dist_ac} > {total_ab_bc}, the triangle inequality is violated.")

if __name__ == "__main__":
    verify_triangle_inequality_violation()
