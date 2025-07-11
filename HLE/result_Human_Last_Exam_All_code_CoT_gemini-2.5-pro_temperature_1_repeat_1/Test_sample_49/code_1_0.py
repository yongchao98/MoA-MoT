def pyramid_f_vector(base_f_vector):
    """
    Calculates the f-vector of a pyramid constructed over a base polytope.
    The base polytope's f-vector is given as a list [f_0, f_1, ...].
    The dimension of the base is len(base_f_vector).
    The pyramid has dimension len(base_f_vector) + 1.
    """
    # The f-vector of a k-polytope has k+1 elements (f_0 to f_k)
    # To handle indexing easily, we treat the polytope itself as f_{dim} = 1
    # and use a conventional f_{-1} = 1 for the new apex.
    # f_i(Pyramid) = f_i(Base) + f_{i-1}(Base)
    
    # Prepend f_{-1}=1 for the apex calculation
    extended_base_f = [1] + base_f_vector
    
    pyramid_f = []
    # f_i of pyramid = f_i of base + f_{i-1} of base
    for i in range(len(extended_base_f)):
        # Get f_i of base. If i is out of bounds, it's 0.
        f_i_base = base_f_vector[i] if i < len(base_f_vector) else 0
        
        # Get f_{i-1} of base from the extended list
        f_im1_base = extended_base_f[i]
        
        pyramid_f.append(f_i_base + f_im1_base)
        
    return pyramid_f

# Step 1: Define the f-vector of a square (a 2-polytope).
# f_0=4 vertices, f_1=4 edges.
square_f_vector = [4, 4]

# Step 2: Calculate the f-vector of the square pyramid (a 3-polytope).
# This is a pyramid over the square.
sq_pyramid_f_vector = pyramid_f_vector(square_f_vector)

# Step 3: Calculate the f-vector of the 4-polytope in question.
# This is a pyramid over the square pyramid.
final_f_vector = pyramid_f_vector(sq_pyramid_f_vector)

# Print the resulting f-vector components as requested.
f0, f1, f2, f3 = final_f_vector

print(f"The f-vector of the non-simplicial 4-polytope with 6 vertices is (f_0, f_1, f_2, f_3).")
print(f"f_0 (vertices) = {f0}")
print(f"f_1 (edges) = {f1}")
print(f"f_2 (2-faces) = {f2}")
print(f"f_3 (3-faces/facets) = {f3}")