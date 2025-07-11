import sys

def count_phosphorus_colors():
    """
    This function identifies and counts the number of distinct colors
    of the main allotropes of phosphorus.
    """
    # The four major allotropes of phosphorus have distinct colors.
    # 1. White Phosphorus (often turns yellow upon exposure to light)
    # 2. Red Phosphorus
    # 3. Violet Phosphorus (Hittorf's phosphorus)
    # 4. Black Phosphorus
    colors = ["White", "Red", "Violet", "Black"]
    
    # Count the number of colors
    total_colors = len(colors)
    
    # As requested, create and print an equation showing how the total was calculated.
    # We will represent each counted color as '1' in the sum.
    equation_parts = ['1' for color in colors]
    equation_str = " + ".join(equation_parts)
    
    print(f"The calculation for the number of colors is:")
    print(f"{equation_str} = {total_colors}")
    
    print(f"\nPhosphorus has {total_colors} major allotropes with distinct colors: {', '.join(colors)}.")

count_phosphorus_colors()