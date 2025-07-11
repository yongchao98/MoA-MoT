import numpy as np

def analyze_flag_rank(name, matrix):
    """Calculates and prints the rank of a given flag matrix."""
    rank = np.linalg.matrix_rank(matrix)
    print(f"--- {name} ---")
    print("Matrix representation:")
    print(matrix)
    print(f"The rank of the flag of {name} is {rank}.\n")
    return rank

def solve_task():
    """
    Identifies African nation flags with the same algebraic rank as Denmark's flag.
    Colors are mapped to distinct integers (e.g., Red=1, White=2, Blue=3, etc.)
    to ensure maximal rank.
    """
    
    # 1. Analyze the flag of Denmark to find the target rank.
    # R=1, W=2
    denmark_matrix = np.array([
        [1, 1, 2, 1, 1],
        [1, 1, 2, 1, 1],
        [2, 2, 2, 2, 2],
        [1, 1, 2, 1, 1],
        [1, 1, 2, 1, 1]
    ])
    target_rank = analyze_flag_rank("Denmark", denmark_matrix)

    print("==============================================")
    print("Analyzing African flags...\n")
    
    # Dictionary to hold flag matrices for African nations
    african_flags = {
        # Rank 2 Candidate: Horizontal bicolor with a central emblem
        "Angola": np.array([
            [1, 1, 3, 1, 1],  # R=1, B=2, Emblem=3
            [1, 1, 3, 1, 1],
            [2, 2, 3, 2, 2],
            [2, 2, 3, 2, 2]
        ]),
        # Rank 2 Candidate: Horizontal tricolor with a central emblem
        "Egypt": np.array([
            [1, 1, 1, 1, 1],  # R=1, W=2, B=3, Eagle=4
            [2, 2, 4, 4, 2],
            [3, 3, 3, 3, 3]
        ]),
        # Rank 2 Candidate: Striped field with a canton
        "Liberia": np.array([
            [1, 1, 2, 2, 2, 2], # BlueCanton=1, WhiteStripe=2, RedStripe=3
            [1, 1, 3, 3, 3, 3],
            [2, 2, 2, 2, 2, 2],
            [3, 3, 3, 3, 3, 3]
        ]),
        # Rank 2 Candidate: Striped field with an emblem in one stripe
        "Rwanda": np.array([
            [1, 1, 4, 1, 1], # Blue=1, Yellow=2, Green=3, Sun=4
            [1, 1, 1, 1, 1],
            [2, 2, 2, 2, 2],
            [3, 3, 3, 3, 3]
        ]),
        # Rank 2 Candidate: Striped field with a canton
        "Togo": np.array([
            [1, 1, 2, 2, 2], # RedCanton=1, GreenStripe=2, YellowStripe=3
            [1, 1, 3, 3, 3],
            [2, 2, 2, 2, 2],
            [3, 3, 3, 3, 3]
        ]),
        # Rank 1 Control Case: Simple vertical stripes
        "Nigeria": np.array([
            [1, 2, 1], # G=1, W=2
            [1, 2, 1],
            [1, 2, 1]
        ]),
        # Rank 1 Control Case: Simple horizontal stripes
        "Gabon": np.array([
            [1, 1, 1, 1], # G=1, Y=2, B=3
            [2, 2, 2, 2],
            [3, 3, 3, 3]
        ])
    }

    matching_flags = []
    for name, matrix in african_flags.items():
        rank = analyze_flag_rank(name, matrix)
        if rank == target_rank:
            matching_flags.append(name)
            
    print("==============================================")
    print(f"The flag of Denmark has a rank of {target_rank}.")
    print("The flags of the following African nations also have a rank of 2:")
    print(", ".join(sorted(matching_flags)))
    
    # Final answer in the specified format
    final_answer = ", ".join(sorted(matching_flags))
    print(f"\n<<<__{final_answer}__>>>")

# Run the analysis
solve_task()