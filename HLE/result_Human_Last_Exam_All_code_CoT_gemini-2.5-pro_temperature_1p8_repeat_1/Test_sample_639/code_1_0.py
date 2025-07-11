import math

def calculate_and_print_distances():
    """
    Calculates and prints the set of possible normalized distances (r)
    between two hard spheres on a plane for r <= 3.
    """
    print("The possible normalized distances r for two hard spheres on a hyperplane are:")

    # r = 1.00: Direct contact
    r1 = 1.0
    print(f"r = {r1:.2f} (Direct contact)")

    # r = sqrt(2) ~ 1.41: Diagonal of a unit square in square packing
    d_sq_2 = 2.0
    r2 = math.sqrt(d_sq_2)
    print(f"r = sqrt({d_sq_2:.2f}) = {r2:.2f} (From square packing, sqrt(1^2+1^2))")

    # r = sqrt(3) ~ 1.73: Second neighbor distance in hexagonal packing
    d_sq_3 = 3.0
    r3 = math.sqrt(d_sq_3)
    print(f"r = sqrt({d_sq_3:.2f}) = {r3:.2f} (From hexagonal packing)")

    # r = 2.00: Collinear contact (i-k-j)
    r4 = 2.0
    print(f"r = {r4:.2f} (Three spheres in a line)")

    # r = sqrt(5) ~ 2.24: From square packing (2,1)
    d_sq_5 = 5.0
    r5 = math.sqrt(d_sq_5)
    print(f"r = sqrt({d_sq_5:.2f}) = {r5:.2f} (From square packing, sqrt(2^2+1^2))")

    # r = sqrt(7) ~ 2.65: From hexagonal packing
    d_sq_7 = 7.0
    r6 = math.sqrt(d_sq_7)
    print(f"r = sqrt({d_sq_7:.2f}) = {r6:.2f} (From hexagonal packing)")

    # r = sqrt(8) ~ 2.83: From square packing (2,2)
    d_sq_8 = 8.0
    r7 = math.sqrt(d_sq_8)
    print(f"r = sqrt({d_sq_8:.2f}) = {r7:.2f} (From square packing, sqrt(2^2+2^2))")

    # r = 3.00: Collinear contact (i-k-l-j)
    r8 = 3.0
    print(f"r = {r8:.2f} (Four spheres in a line)")

calculate_and_print_distances()
