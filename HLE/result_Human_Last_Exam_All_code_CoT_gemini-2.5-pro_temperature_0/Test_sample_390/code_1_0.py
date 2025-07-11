import numpy as np

def calculate_s_point(s, y_vectors):
    """Calculates a point in the set S for a given s."""
    s = np.array(s)
    point = []
    for y in y_vectors:
        y = np.array(y)
        inner_product = np.dot(y, s)
        point.append(inner_product**2)
    return tuple(point)

def main():
    """
    Demonstrates that for n=3, the set S is not a simplex (i.e., not planar).
    """
    # Define the linearly independent vectors for n=3
    y1 = [1, 0, 0]
    y2 = [0, 1, 0]
    y3 = [1, 1, 1]
    y_vectors = [y1, y2, y3]

    # Choose some unit vectors s
    s1 = [1, 0, 0]
    s2 = [0, 1, 0]
    s3 = [0, 0, 1]
    s4 = [1/np.sqrt(2), 1/np.sqrt(2), 0]

    # Calculate the corresponding points in S
    p1 = calculate_s_point(s1, y_vectors)
    p2 = calculate_s_point(s2, y_vectors)
    p3 = calculate_s_point(s3, y_vectors)
    p4 = calculate_s_point(s4, y_vectors)

    print(f"For s = {s1}, the point in S is ({p1[0]:.2f}, {p1[1]:.2f}, {p1[2]:.2f})")
    print(f"For s = {s2}, the point in S is ({p2[0]:.2f}, {p2[1]:.2f}, {p2[2]:.2f})")
    print(f"For s = {s3}, the point in S is ({p3[0]:.2f}, {p3[1]:.2f}, {p3[2]:.2f})")
    print("These three points lie on the plane x_3 = 1.")
    
    print(f"\nNow, for s = ({s4[0]:.2f}, {s4[1]:.2f}, {s4[2]:.2f}), the point in S is ({p4[0]:.2f}, {p4[1]:.2f}, {p4[2]:.2f})")
    print("This point's third component is not 1, so it does not lie on the same plane.")
    print("This shows the shape is not a simplex.")

if __name__ == "__main__":
    main()