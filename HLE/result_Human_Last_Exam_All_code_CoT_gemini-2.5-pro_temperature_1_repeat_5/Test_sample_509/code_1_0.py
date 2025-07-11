import numpy as np

def section_map_for_open_square(p1):
    """
    Constructs a section for the fibration pi_{1,2} on the open unit square M = (0,1)^2.
    
    The map is pi_{1,2}: conf_2(M) -> conf_1(M).
    A section s: conf_1(M) -> conf_2(M) must satisfy pi_{1,2}(s(p)) = p.
    This means s(p1) must be of the form (p1, p2), where p2 is a point in M
    that depends continuously on p1 and is never equal to p1.

    Args:
        p1 (tuple): A point (x1, y1) in the configuration space conf_1(M),
                    which is just M itself. We assume 0 < x1 < 1 and 0 < y1 < 1.

    Returns:
        tuple: A configuration (x1, y1, x2, y2) in conf_2(M).
    """
    x1, y1 = p1
    
    # Check if the input point is in the open unit square
    if not (0 < x1 < 1 and 0 < y1 < 1):
        raise ValueError("Input point must be in the open unit square (0,1)x(0,1)")

    # Define the second point p2 = (x2, y2) using a continuous function of p1.
    # We choose a simple function that guarantees p2 is in M and p2 != p1.
    # Let x2 = (x1 + 1) / 2 and y2 = (y1 + 1) / 2.
    # For any x1 in (0,1), x2 will be in (0.5, 1), so x2 > x1.
    # For any y1 in (0,1), y2 will be in (0.5, 1), so y2 > y1.
    # Therefore, p2 is always different from p1 and always inside M.
    x2 = (x1 + 1) / 2
    y2 = (y1 + 1) / 2
    
    p2 = (x2, y2)
    
    # The resulting configuration in conf_2(M)
    config_2 = (x1, y1, x2, y2)
    
    return config_2

def main():
    """
    Main function to demonstrate the section map.
    """
    print("Demonstrating the construction of a section for pi_{1,2} on M = (0,1)^2.")
    print("The section maps a point p1=(x1, y1) to a configuration (p1, p2).")
    print("We use the equation p2 = ((x1+1)/2, (y1+1)/2) to define the second point.")
    
    # Example point
    p1_example = (0.2, 0.3)
    
    # Calculate the section
    result_config = section_map_for_open_square(p1_example)
    
    x1, y1, x2, y2 = result_config
    
    print("\n--- Example ---")
    print(f"Input point p1 = ({x1}, {y1})")
    
    print("\nEquation for the second point p2=(x2, y2):")
    print(f"x2 = (x1 + 1) / 2 = ({x1} + 1) / 2 = {x2}")
    print(f"y2 = (y1 + 1) / 2 = ({y1} + 1) / 2 = {y2}")
    
    print(f"\nResulting point p2 = ({x2}, {y2})")
    
    print(f"\nFinal configuration in conf_2(M) is ({x1}, {y1}, {x2}, {y2})")
    
    # Verify the properties
    is_distinct = (x1, y1) != (x2, y2)
    is_p2_in_M = (0 < x2 < 1) and (0 < y2 < 1)
    
    print("\nVerification:")
    print(f"Is p1 distinct from p2? {is_distinct}")
    print(f"Is p2 inside M=(0,1)^2? {is_p2_in_M}")


if __name__ == "__main__":
    main()

<<<E>>>