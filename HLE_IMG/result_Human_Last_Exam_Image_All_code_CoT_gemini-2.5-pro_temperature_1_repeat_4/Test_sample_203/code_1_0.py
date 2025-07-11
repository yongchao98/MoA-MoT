def find_heesch_numbers():
    """
    This function determines the Heesch number for the three identical polygons
    shown in the image.

    The reasoning is as follows:
    1.  The polygon's shape is centrally symmetric. This can be observed by noting
        that for every feature (like a protrusion or indentation) on its boundary,
        an identical feature exists on the opposite side, rotated by 180 degrees.

    2.  A fundamental theorem of geometry states that any centrally symmetric polygon
        can tile the plane. This means copies of the shape can cover a flat surface
        without gaps or overlaps.

    3.  The Heesch number of a shape that can tile the plane is defined to be
        infinity. This is because it can be surrounded by an endless number of layers
        of its copies.

    4.  The image shows three identical polygons. Therefore, the Heesch number for
        each is infinity.
    """
    
    # The infinity symbol (∞) is the answer for a tile that tessellates the plane.
    heesch_number_for_one_polygon = u'\u221e'
    
    # There are three identical polygons in the figure.
    number_of_polygons = 3
    
    # Create a list of the answers.
    results = [heesch_number_for_one_polygon] * number_of_polygons
    
    # Join the results with commas for the final output.
    final_answer = ", ".join(results)
    
    print("The Heesch number for a polygon that tiles the plane is infinity.")
    print(f"The polygon in the image is centrally symmetric, so it tiles the plane.")
    print(f"Therefore, its Heesch number is {heesch_number_for_one_polygon}.")
    print("\nSince there are three identical polygons, the Heesch numbers are:")
    print(final_answer)

find_heesch_numbers()