import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

def compute_rank(X, svd_rank=0):
    U, s, _ = np.linalg.svd(X, full_matrices=False)

    def omega(x):
        return 0.56 * x**3 - 0.95 * x**2 + 1.82 * x + 1.43

    if svd_rank == 0:
        beta = np.divide(*sorted(X.shape))
        tau = np.median(s) * omega(beta)
        rank = np.sum(s > tau)
        if rank == 0:
            rank = 1
    elif 0 < svd_rank < 1:
        cumulative_energy = np.cumsum(s**2 / (s**2).sum())
        rank = np.searchsorted(cumulative_energy, svd_rank) + 1
    elif svd_rank >= 1 and isinstance(svd_rank, int):
        rank = min(svd_rank, U.shape[1])
    else:
        rank = min(X.shape)

    return rank

def compute_svd(X, svd_rank=0):
    rank = compute_rank(X, svd_rank)
    U, s, V = np.linalg.svd(X, full_matrices=False)
    V = V.conj().T

    U = U[:, :rank]
    V = V[:, :rank]
    s = s[:rank]

    return U, s, V

def pseudo_hankel_matrix(X: np.ndarray, d: int):
    return (
        sliding_window_view(X.T, (d, X.shape[0]))[:, 0]
        .reshape(X.shape[1] - d + 1, -1)
        .T
    )

def main_solution(matrix_data, svd_rank, hankel_depth):
    X = np.array(matrix_data)
    hankel_matrix = pseudo_hankel_matrix(X, hankel_depth)
    U, _, _ = compute_svd(hankel_matrix, svd_rank)
    return U.tolist()

# Given input
matrix_data = [[-0.9991242666198712, -8.84472174920133, -5.446824505595145, -0.03918737823216567, -4.024622505503695, 8.33109386546504, -8.437119604657523, -5.291314206871169, 0.11970419812017141, -1.4393642784704639], [3.3064976583552337, 8.247323019691105, -3.8362845952013913, 4.563910148753795, 7.582710313255699, -3.6661110720550516, 0.4410700758225339, 2.570761706301436, -9.824449758419204, -3.256125059373609], [-5.944379696650408, 1.9691164380588155, -6.401640923428893, -1.8970287036935591, 3.2578691341977972, 7.247219450113974, -7.114236440701673, -2.7387752900515983, -6.314198368656065, 3.7655505490043737], [-9.471865535266266, 7.648283897198301, 2.7199437445493118, -7.57397473108675, 6.834724085653903, -4.675077246147923, -5.669398787020663, -0.8799422493423066, 1.5021314693549712, 6.366379528619468], [-7.733389487171111, 7.626931227592955, -4.665639033122389, -2.958831401752276, 5.288176338656507, -0.3825263781010353, -9.410711410255184, -8.794620421758658, -8.485279309903795, -3.832395854522961], [9.94767781447036, -6.238617694241797, 6.293882041662844, -5.450366346287405, -4.901723788001782, -9.64671333606017, -0.1269090392170984, 0.27180052066580984, -3.7790008517075435, 8.389977246555013]]
svd_rank = 5
hankel_depth = 10

# Compute the result
result = main_solution(matrix_data, svd_rank, hankel_depth)
print(result)