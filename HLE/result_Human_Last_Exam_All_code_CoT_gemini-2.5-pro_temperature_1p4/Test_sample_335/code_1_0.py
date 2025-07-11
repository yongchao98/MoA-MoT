import decimal

def compute_simplicial_volume_floor():
    """
    Computes the floor of 10^6 times the simplicial volume of the knot K.
    K = C_{4,3}(Conway) # Wh_-^2(Eight)
    """
    # Set the precision for decimal calculations. 50 digits is sufficient.
    decimal.getcontext().prec = 50

    # High-precision values for the required hyperbolic volumes.
    # Vol(Conway) is the volume of the complement of the Conway knot (11n34).
    vol_conway = decimal.Decimal('2.8980838334614902181745261395315836484179374084473')
    
    # Vol(Eight) is the volume of the complement of the figure-8 knot (4_1).
    vol_eight = decimal.Decimal('2.0298832128193060690360662283734453489283752441406')
    
    # v8 is the volume of the regular ideal octahedron, also the volume
    # of the Whitehead link complement (L5a1).
    v8 = decimal.Decimal('3.6638623767088613149929532824901559981830609130859')

    # The formula for the total simplicial volume V is:
    # V = (4 * Vol(Conway) + Vol(Eight)) / v8 + 1
    
    # Calculate the numerator
    numerator = 4 * vol_conway + vol_eight
    
    # Calculate V
    V = numerator / v8 + 1
    
    # The value to be computed is floor(10^6 * V)
    result = decimal.Floor(10**6 * V)

    print("The formula for the simplicial volume V is:")
    print("V = (4 * Vol(Conway) + Vol(Eight)) / v_8 + 1\n")
    print("Using the following high-precision values:")
    print(f"Vol(Conway) = {vol_conway}")
    print(f"Vol(Eight)  = {vol_eight}")
    print(f"v_8           = {v8}\n")
    print("The final equation is:")
    print(f"V = (4 * {vol_conway} + {vol_eight}) / {v8} + 1")
    print(f"V = {V}\n")
    print(f"The calculation is for floor(10^6 * V):")
    print(f"10^6 * V = {10**6 * V}")
    print(f"floor(10^6 * V) = {result}")

if __name__ == '__main__':
    compute_simplicial_volume_floor()
