import numpy as np

def create_tonga_flag_matrix(height, width, a, b):
    """
    Creates a matrix representation of the flag of Tonga.
    
    The flag is `height` x `width`.
    Red pixels have value `a`, white pixels have value `b`.
    """
    # Initialize the entire flag as red (value 'a')
    flag = np.full((height, width), a, dtype=float)
    
    # Define the canton (white rectangle)
    # Let's say canton is 1/2 height and 1/3 width
    canton_h = height // 2
    canton_w = width // 3
    flag[0:canton_h, 0:canton_w] = b
    
    # Define the red cross within the canton
    # Make the cross proportional to the canton size
    cross_h_start = canton_h // 4
    cross_h_end = cross_h_start * 3
    cross_w_start = canton_w // 4
    cross_w_end = cross_w_start * 2 # a bit thinner vertically

    # Horizontal bar of the cross
    flag[cross_h_start:cross_h_end, 0:canton_w] = a
    # Vertical bar of the cross
    flag[0:canton_h, cross_w_start:cross_w_end] = a
    
    return flag

def main():
    """
    Calculates and prints the maximal rank of the Tonga flag matrix.
    """
    # To find the maximal rank, we choose a and b to be distinct non-zero numbers.
    a = 1.0
    b = 2.0
    
    # Define flag dimensions (must be large enough to have all features)
    height = 60
    width = 90
    
    # Create the matrix
    tonga_matrix = create_tonga_flag_matrix(height, width, a, b)
    
    # Calculate the rank
    rank = np.linalg.matrix_rank(tonga_matrix)
    
    print(f"Assigning red pixels the value a={a} and white pixels the value b={b}.")
    print("The structure of the flag of Tonga results in a matrix with two types of unique rows.")
    print("1. Rows that are all red (value 'a').")
    print("2. Rows that pass through the white canton and the vertical bar of the cross (values 'a' and 'b').")
    print("These two row types are linearly independent when a != b.")
    print(f"The maximal possible rank of this matrix is therefore 2.")
    print(f"\nVerification with a sample {height}x{width} matrix:")
    print(f"The calculated rank is: {rank}")

if __name__ == "__main__":
    main()
<<<2>>>