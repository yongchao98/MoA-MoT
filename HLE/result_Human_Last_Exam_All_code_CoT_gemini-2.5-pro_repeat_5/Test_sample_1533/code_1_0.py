def calculate_bm_mi_ratio(a, b, c):
    """
    Calculates the ratio BM/MI for a triangle with side lengths a, b, c.
    
    Args:
        a (float): The length of side BC.
        b (float): The length of side AC.
        c (float): The length of side AB.
        
    Returns:
        float: The value of the ratio BM/MI.
    """
    if a + b <= c or a + c <= b or b + c <= a:
        print("Error: The given side lengths do not form a valid triangle.")
        return None
        
    # The derived formula for the ratio BM/MI is (a+c)/b
    ratio = (a + c) / b
    return ratio

# --- Example Usage ---
# Let's use a sample triangle with side lengths a=7, b=8, c=9.
a = 7
b = 8
c = 9

# Calculate the ratio
result = calculate_bm_mi_ratio(a, b, c)

if result is not None:
    # Print the final result, showing the formula and the numbers used.
    print(f"For a triangle with side lengths a = {a}, b = {b}, c = {c}:")
    print("The ratio BM/MI is expressed by the formula: (a + c) / b")
    print(f"Substituting the values: ({a} + {c}) / {b}")
    print(f"= {a+c} / {b}")
    print(f"= {result}")
