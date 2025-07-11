import numpy as np

def demonstrate_section():
    """
    Illustrates the concept of a section for a projection between
    configuration spaces on a non-compact manifold (M = R).
    """
    # A point in conf_1(R) is just a single number, represented as a list.
    point_x1 = [17.0]

    # A section s: conf_1(R) -> conf_2(R) must continuously map a point x1
    # to a distinct pair (x1, x2). We can define a simple continuous rule.
    # For a given point p = [x1], the section returns [x1, x1 + 1].
    section_map = lambda p: [p[0], p[0] + 1]

    # The projection pi_1,2: conf_2(R) -> conf_1(R) forgets the second point.
    projection_map = lambda pair: [pair[0]]

    print(f"Let's start with a point in conf_1(R): {point_x1}")

    # Apply the section map s
    point_in_conf2 = section_map(point_x1)
    print(f"Applying the section s gives a point in conf_2(R): {point_in_conf2}")

    # Apply the projection map pi_1,2 to the result of the section
    result = projection_map(point_in_conf2)
    print(f"Applying the projection pi_1,2 to that result gives: {result}")

    # For a strict section, the result must be the original point.
    is_identity = (result == point_x1)
    print(f"This demonstrates that (pi_1,2 o s) is the identity map: {is_identity}")

demonstrate_section()
