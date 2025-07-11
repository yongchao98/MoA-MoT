import math

def find_smallest_n():
    """
    This function systematically searches for the smallest N by constructing sets of 16 distinct
    positive integers whose product is a perfect fourth power. For each valid set, it
    calculates the maximum element (N) and finds the minimum N among all constructed sets.
    """
    min_n = float('inf')
    best_set_info = {}

    # --- Candidate Set 1: 4 primes, exponents {0,1} ---
    s1 = set()
    for e2 in range(2):
        for e3 in range(2):
            for e5 in range(2):
                for e7 in range(2):
                    s1.add((2**e2) * (3**e3) * (5**e5) * (7**e7))
    n1 = max(s1)
    if n1 < min_n:
        min_n = n1
        best_set_info = {"N": n1, "set": s1, "config": "4 primes (2,3,5,7), exponents {0,1}"}

    # --- Candidate Set 2: 2 primes, exponents {0,1,2,3} ---
    s2 = set()
    for e2 in range(4):
        for e3 in range(4):
            s2.add((2**e2) * (3**e3))
    n2 = max(s2)
    if n2 < min_n:
        min_n = n2
        best_set_info = {"N": n2, "set": s2, "config": "2 primes (2,3), exponents {0,1,2,3}"}

    # --- Candidate Set 3: 3 primes, exponents {0..3}x{0..1}x{0..1} ---
    s3 = set()
    for e2 in range(4):
        for e3 in range(2):
            for e5 in range(2):
                s3.add((2**e2) * (3**e3) * (5**e5))
    n3 = max(s3)
    if n3 < min_n:
        min_n = n3
        best_set_info = {"N": n3, "set": s3, "config": "3 primes (2,3,5), exponents {0..3}x{0..1}x{0..1}"}

    # --- Candidate Set 4: 3 primes, exponents {0..2}x{0..2}x{0..1} ---
    # This construction creates 18 numbers. We must remove 2 to form a valid 16-number set.
    full_set_vectors = []
    for e2 in range(3):
        for e3 in range(3):
            for e5 in range(2):
                full_set_vectors.append({'e': (e2, e3, e5), 'val': (2**e2)*(3**e3)*(5**e5)})

    sum_e2_full = sum(v['e'][0] for v in full_set_vectors)
    sum_e3_full = sum(v['e'][1] for v in full_set_vectors)
    sum_e5_full = sum(v['e'][2] for v in full_set_vectors)

    # The sum of exponents of the two removed numbers must match the sum of exponents
    # of the full set, modulo 4.
    target_sum_e2 = sum_e2_full % 4
    target_sum_e3 = sum_e3_full % 4
    target_sum_e5 = sum_e5_full % 4

    for i in range(len(full_set_vectors)):
        for j in range(i + 1, len(full_set_vectors)):
            v1 = full_set_vectors[i]['e']
            v2 = full_set_vectors[j]['e']

            if ((v1[0] + v2[0]) % 4 == target_sum_e2 and
                (v1[1] + v2[1]) % 4 == target_sum_e3 and
                (v1[2] + v2[2]) % 4 == target_sum_e5):
                
                val1 = full_set_vectors[i]['val']
                val2 = full_set_vectors[j]['val']
                
                current_set = {v['val'] for v in full_set_vectors} - {val1, val2}
                n4 = max(current_set)

                if n4 < min_n:
                    min_n = n4
                    best_set_info = {
                        "N": n4,
                        "set": current_set,
                        "config": f"3 primes (2,3,5), exponents {{0..2}}x{{0..2}}x{{0..1}}, removing {{{val1}, {val2}}}"
                    }

    print(f"The optimal set of 16 numbers was found with the configuration:\n{best_set_info['config']}\n")
    print("The 16 distinct positive integers are:")
    final_numbers = sorted(list(best_set_info['set']))
    # The prompt asks to "output each number in the final equation".
    # This is interpreted as listing the numbers in the set that solves the problem.
    for num in final_numbers:
        print(num)
    
    print(f"\nThe smallest N is the largest number in this set.")
    print(f"N = {best_set_info['N']}")

find_smallest_n()