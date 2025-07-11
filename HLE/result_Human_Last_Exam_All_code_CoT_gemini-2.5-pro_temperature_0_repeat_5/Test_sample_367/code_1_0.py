import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def solve():
    """
    Finds the number of distinct parallelograms satisfying the given conditions.
    """
    # 1. Generate Pythagorean vectors (c,d) where hypotenuse e < 100.
    py_vectors = set()
    limit = 100
    # Use Euclid's formula to generate Pythagorean triples (c,d,e)
    for m in range(2, int(math.sqrt(limit)) + 1):
        for n in range(1, m):
            if (m - n) % 2 == 1 and gcd(m, n) == 1:
                # Primitive triple
                c_prim = m*m - n*n
                d_prim = 2*m*n
                e_prim = m*m + n*n
                
                # Generate all non-primitive triples by scaling with k
                k = 1
                while k * e_prim < limit:
                    c, d, e = k * c_prim, k * d_prim, k * e_prim
                    # Add all orientations of the vector
                    py_vectors.add((c, d))
                    py_vectors.add((d, c))
                    py_vectors.add((-c, d))
                    py_vectors.add((d, -c))
                    py_vectors.add((c, -d))
                    py_vectors.add((-d, c))
                    py_vectors.add((-c, -d))
                    py_vectors.add((-d, -c))
                    k += 1
    
    vector_list = list(py_vectors)
    found_parallelograms = set()

    # 2. Iterate through all unique pairs of vectors
    for i in range(len(vector_list)):
        for j in range(i + 1, len(vector_list)):
            v1 = vector_list[i]
            v2 = vector_list[j]
            
            x, y = v1
            z, w = v2
            
            # 3a. Calculate side lengths
            a_sq = x*x + y*y
            b_sq = z*z + w*w
            a = int(math.sqrt(a_sq))
            b = int(math.sqrt(b_sq))

            # 3b. Check side length conditions
            if a == b:  # Not a rhombus
                continue
            if a + b >= 100:
                continue
            if gcd(a, b) != 1:  # Coprime
                continue

            # 3c. Calculate diagonal lengths
            d1_sq = (x + z)**2 + (y + w)**2
            d2_sq = (x - z)**2 + (y - w)**2
            
            # The vectors are already Pythagorean, so their sums/diffs must also be.
            # This check is implicitly handled by the construction method.
            # We just need the integer values for the key.
            if not is_perfect_square(d1_sq) or not is_perfect_square(d2_sq):
                continue
            d1 = int(math.sqrt(d1_sq))
            d2 = int(math.sqrt(d2_sq))

            # 3d. Check non-rectangle condition
            if x*z + y*w == 0:
                continue

            # 4. Store canonical representation of the parallelogram
            sides = tuple(sorted((a, b)))
            diagonals = tuple(sorted((d1, d2)))
            key = (sides, diagonals)
            found_parallelograms.add(key)

    # 5. Output the result
    print(f"Found {len(found_parallelograms)} distinct parallelograms.")

solve()
<<<8>>>