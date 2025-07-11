import math

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __sub__(self, other):
        return Point(self.x - other.x, self.y - other.y)

    def __repr__(self):
        return f"({self.x:.4f}, {self.y:.4f})"

def solve_square():
    """
    Finds the vertices of a square given one point on each side.
    """
    points = [
        Point(0.3511, 0.2027),
        Point(0.6753, 0.8303),
        Point(-0.2845, 0.9905),
        Point(-0.128, 0.2218)
    ]
    
    # There are 3 ways to pair the 4 points as being on opposite sides.
    # The order of points in the tuple defines the consecutive order on the square's perimeter.
    # e.g., in (0, 1, 2, 3), P0 is on side V0V1, P1 on V1V2, P2 on V2V3, P3 on V3V0
    point_orderings = [
        (0, 1, 2, 3), # P0-P2 opposite, P1-P3 opposite
        (0, 2, 1, 3), # P0-P1 opposite, P2-P3 opposite
        (0, 1, 3, 2)  # P0-P3 opposite, P1-P2 opposite
    ]

    for order in point_orderings:
        p1, p2, p3, p4 = points[order[0]], points[order[1]], points[order[2]], points[order[3]]

        # p1, p3 are on opposite sides; p2, p4 are on opposite sides
        vec1 = p1 - p3
        vec2 = p2 - p4
        
        u, v = vec1.x, vec1.y
        p, q = vec2.x, vec2.y
        
        # Calculate the two possible slopes for the sides
        slopes = []
        if abs(u - q) > 1e-9:
            slopes.append((v + p) / (u - q))
        if abs(u + q) > 1e-9:
            slopes.append((v - p) / (u + q))
            
        for m in slopes:
            if abs(m) < 1e-9: continue
            m_perp = -1.0 / m
            
            # Lines L1 (thru p1), L3 (thru p3) have slope m
            # Lines L2 (thru p2), L4 (thru p4) have slope m_perp
            # Line equation y = slope * x + d
            d1 = p1.y - m * p1.x
            d3 = p3.y - m * p3.x
            d2 = p2.y - m_perp * p2.x
            d4 = p4.y - m_perp * p4.x

            m_diff = m - m_perp
            if abs(m_diff) < 1e-9: continue

            # Vertices are intersections of the lines
            # V1 = L1 intersect L4
            x_v1 = (d4 - d1) / m_diff
            y_v1 = m * x_v1 + d1
            v1 = Point(x_v1, y_v1)

            # V2 = L1 intersect L2
            x_v2 = (d2 - d1) / m_diff
            y_v2 = m * x_v2 + d1
            v2 = Point(x_v2, y_v2)

            # V3 = L3 intersect L2
            x_v3 = (d2 - d3) / m_diff
            y_v3 = m_perp * x_v3 + d2
            v3 = Point(x_v3, y_v3)

            # V4 = L3 intersect L4
            x_v4 = (d4 - d3) / m_diff
            y_v4 = m * x_v4 + d3
            v4 = Point(x_v4, y_v4)
            
            # Check if points lie on the segments
            tol = 1e-6
            p1_on_v1v2 = (min(v1.x, v2.x) - tol <= p1.x <= max(v1.x, v2.x) + tol and
                          min(v1.y, v2.y) - tol <= p1.y <= max(v1.y, v2.y) + tol)
            p2_on_v2v3 = (min(v2.x, v3.x) - tol <= p2.x <= max(v2.x, v3.x) + tol and
                          min(v2.y, v3.y) - tol <= p2.y <= max(v2.y, v3.y) + tol)
            p3_on_v3v4 = (min(v3.x, v4.x) - tol <= p3.x <= max(v3.x, v4.x) + tol and
                          min(v3.y, v4.y) - tol <= p3.y <= max(v3.y, v4.y) + tol)
            p4_on_v4v1 = (min(v4.x, v1.x) - tol <= p4.x <= max(v4.x, v1.x) + tol and
                          min(v4.y, v1.y) - tol <= p4.y <= max(v4.y, v1.y) + tol)

            if p1_on_v1v2 and p2_on_v2v3 and p3_on_v3v4 and p4_on_v4v1:
                vertices = sorted([v1, v2, v3, v4], key=lambda pt: pt.x)
                result_str = []
                for v in vertices:
                    # Print each coordinate pair of the solution
                    print(f"({v.x:.2f}, {v.y:.2f})")
                    result_str.append(f"({v.x:.2f}, {v.y:.2f})")
                
                # Final answer format
                # print(f"<<<{' '.join(result_str)}>>>")
                return

if __name__ == '__main__':
    solve_square()
