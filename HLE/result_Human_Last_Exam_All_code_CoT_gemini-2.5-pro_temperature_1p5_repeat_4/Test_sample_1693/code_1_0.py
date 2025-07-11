def find_largest_sunk_ship():
    """
    This function provides information on the largest French ship by displacement 
    sunk by a U-boat before the 1940 armistice.
    The information is based on historical research.
    """
    
    # Ship details
    ship_name = "Champlain"
    
    # Displacement is the measure of weight, which was requested.
    # For a large passenger liner, this figure is often close to its Gross Register Tonnage (GRT).
    # Its displacement was 27,680 tons. Its GRT was 28,124.
    displacement = 27680
    
    # Details of the sinking
    cause = "a mine laid by the German submarine U-66"
    date = "June 17, 1940"

    # Print the final answer, ensuring all key details and numbers are included.
    print(f"The largest French ship by displacement sunk by a U-boat before the 1940 armistice was the passenger liner '{ship_name}'.")
    print(f"It was sunk on {date} by {cause}.")
    print(f"The displacement of the '{ship_name}' was {displacement} tons.")

if __name__ == '__main__':
    find_largest_sunk_ship()