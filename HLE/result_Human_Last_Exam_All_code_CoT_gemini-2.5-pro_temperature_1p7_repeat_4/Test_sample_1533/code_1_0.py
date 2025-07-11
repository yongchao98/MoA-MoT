def calculate_bm_mi_ratio(a, b, c):
    """
    Calculates the ratio BM/MI in a triangle ABC.

    In triangle ABC, I is the incenter and M is the intersection of the angle
    bisector BI with the circumcircle. The ratio of the lengths of the segments
    BM to MI can be expressed in terms of the side lengths a, b, and c.

    The derivation shows that BM/MI = (a + c) / b.

    Args:
        a (float): The length of side BC.
        b (float): The length of side CA.
        c (float): The length of side AB.
    
    Returns:
        float: The calculated ratio BM/MI, or None if the side lengths
               do not form a valid triangle.
    """
    # Check if the side lengths form a valid triangle (Triangle Inequality)
    if a + b <= c or a + c <= b or b + c <= a:
        print("Error: The given side lengths do not form a valid triangle.")
        return None

    # The ratio is given by the simple formula (a + c) / b
    ratio = (a + c) / b
    
    # Print the equation with the numbers plugged in
    print(f"The ratio BM/MI is expressed as (a + c) / b.")
    print(f"For the given side lengths:")
    print(f"({a} + {c}) / {b} = {ratio}")
    
    return ratio

# --- Example Usage ---
# You can replace these values with the side lengths of any triangle.

# Example 1: A 3-4-5 right triangle. Let a=3, c=4, so b=5.
print("--- Example 1: a=3, b=5, c=4 ---")
calculate_bm_mi_ratio(a=3, b=5, c=4)

print("\n" + "="*30 + "\n")

# Example 2: An equilateral triangle. a=b=c=6.
print("--- Example 2: a=6, b=6, c=6 ---")
calculate_bm_mi_ratio(a=6, b=6, c=6)

print("\n" + "="*30 + "\n")

# Example 3: A general triangle. a=7, b=8, c=9.
print("--- Example 3: a=7, b=8, c=9 ---")
calculate_bm_mi_ratio(a=7, b=8, c=9)
