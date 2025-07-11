def solve_icosahedron_problem():
    """
    Determines and explains the shape of the water surface in a half-filled
    icosahedron tank standing on one of its faces.
    """
    print("Step 1: Understanding the setup.")
    print("An icosahedron is a polyhedron with 20 identical equilateral triangular faces.")
    print("The tank is standing on one face, so its base is a horizontal triangle.")
    print("-" * 20)

    print("Step 2: Applying symmetry.")
    print("Due to the icosahedron's symmetry, there is a parallel triangular face at the top.")
    print("When the tank is half-filled, the water level is at the horizontal plane exactly halfway between the base and the top.")
    print("This mid-plane is a plane of symmetry for the icosahedron.")
    print("-" * 20)

    print("Step 3: Determining the shape of the water surface.")
    print("The shape of the surface is the cross-section of the icosahedron cut by this horizontal mid-plane.")
    print("This plane intersects the faces forming the 'waist' of the icosahedron.")
    print("The intersection of the plane and the icosahedron forms a polygon.")
    print("-" * 20)

    print("Step 4: Identifying the final shape.")
    print("The mid-plane cuts through 6 of the icosahedron's faces.")
    print("The resulting cross-section is a 6-sided polygon (a hexagon).")
    print("Because of the icosahedron's high degree of symmetry around the vertical axis, this shape is not just any hexagon, but a regular hexagon.")
    print("-" * 20)

    print("Final Answer: The shape of the water surface is a regular hexagon.")

if __name__ == "__main__":
    solve_icosahedron_problem()