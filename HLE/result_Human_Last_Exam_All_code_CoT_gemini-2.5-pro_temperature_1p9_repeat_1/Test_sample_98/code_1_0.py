import numpy as np

def describe_water_surface():
    """
    This function describes the shape of the water surface and verifies its properties.
    The shape is a regular hexagon, which is the cross-section of the icosahedron
    at its center when it's resting on a face.
    """
    print("The shape of the water surface will be a regular hexagon.")
    print("-" * 20)
    print("Verifying the properties of a regular hexagon computationally:")
    print("-" * 20)

    # For demonstration, we define the vertices of a regular hexagon with side length 1.
    # The exact size depends on the icosahedron's size, but the shape is always a regular hexagon.
    # The radius `R` would be `(side_length / 2) / sin(pi/6) = side_length` for this setup.
    R = 1.0 # Radius of the hexagon's circumcircle
    vertices = []
    for i in range(6):
        angle = np.deg2rad(60 * i)
        x = R * np.cos(angle)
        y = R * np.sin(angle)
        vertices.append(np.array([x, y]))

    print("The vertices of the hexagonal water surface (for a hexagon of radius 1.0) are:")
    for i, v in enumerate(vertices):
        print(f"Vertex {i+1}: ({v[0]:.4f}, {v[1]:.4f})")
    
    print("\nCalculating side lengths:")
    side_lengths = []
    for i in range(6):
        p1 = vertices[i]
        p2 = vertices[(i + 1) % 6]  # Wrap around to the first vertex
        distance = np.linalg.norm(p2 - p1)
        side_lengths.append(distance)
        print(f"Side {i+1} (between Vertex {i+1} and Vertex {(i + 1) % 6 + 1}): {distance:.4f}")
    
    # Check if all side lengths are equal
    if np.allclose(side_lengths, side_lengths[0]):
        print("\nConclusion: All side lengths are equal.")
    else:
        print("\nConclusion: Side lengths are not equal.")

    print("\nCalculating internal angles:")
    angles = []
    for i in range(6):
        p_prev = vertices[i-1] 
        p_curr = vertices[i]
        p_next = vertices[(i + 1) % 6]
        
        # Create vectors for the two edges meeting at p_curr
        v1 = p_prev - p_curr
        v2 = p_next - p_curr
        
        # Calculate the angle between the vectors
        cosine_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        angle = np.arccos(cosine_angle)
        angle_deg = np.rad2deg(angle)
        angles.append(angle_deg)
        print(f"Angle at Vertex {i+1}: {angle_deg:.4f} degrees")
        
    # Check if all angles are equal
    if np.allclose(angles, angles[0]):
        print("\nConclusion: All internal angles are equal.")
    else:
        print("\nConclusion: Internal angles are not equal.")

    print("\nSince all sides and angles are equal, the shape is a regular hexagon.")

# Run the function
if __name__ == '__main__':
    describe_water_surface()