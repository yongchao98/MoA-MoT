import math

def triangle_area(x1, y1, z1, x2, y2, z2, x3, y3, z3):
    a = x2*y3 + x3*y1 + x1*y2 - x3*y2 - x1*y3 - x2*y1
    b = x2*z3 + x3*z1 + x1*z2 - x3*z2 - x1*z3 - x2*z1
    c = y2*z3 + y3*z1 + y1*z2 - y3*z2 - y1*z3 - y2*z1
    return 0.5 * math.sqrt(a**2 + b**2 + c**2)

# Coordinates of the points
C = (6, 6*math.sqrt(3), 0)
X = (3, 3*math.sqrt(3), 0)
Y = (9, 3*math.sqrt(3), 0)
Z = (6, 6*math.sqrt(3), 8)

# Calculate the area of each face
area_CXY = triangle_area(*C, *X, *Y)
area_CYZ = triangle_area(*C, *Y, *Z)
area_CXZ = triangle_area(*C, *X, *Z)
area_XYZ = triangle_area(*X, *Y, *Z)

# Total surface area of the tetrahedron CXYZ
total_surface_area = area_CXY + area_CYZ + area_CXZ + area_XYZ

print(total_surface_area)