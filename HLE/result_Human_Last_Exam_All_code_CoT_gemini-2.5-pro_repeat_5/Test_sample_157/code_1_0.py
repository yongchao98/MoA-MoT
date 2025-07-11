import math
import base64

def cosine_similarity(v1, v2):
    """Computes the cosine similarity between two vectors."""
    dot_product = sum(x * y for x, y in zip(v1, v2))
    mag_v1 = math.sqrt(sum(x * x for x in v1))
    mag_v2 = math.sqrt(sum(x * x for x in v2))
    if mag_v1 == 0 or mag_v2 == 0:
        return 0.0
    return dot_product / (mag_v1 * mag_v2)

def main():
    # Define the item vectors
    items = {
        'A': (1, 0, 1), 'B': (0, 1, 1), 'C': (1, 1, 0), 'D': (0, 0, 1),
        'E': (1, 1, 1), 'F': (0, 0, 0), 'G': (1, 0, 0), 'H': (0, 1, 0),
        'I': (1, 0, 1), 'J': (0, 1, 1), 'K': (1, 1, 0), 'L': (0, 0, 1),
        'M': (1, 1, 1), 'N': (0, 0, 0), 'O': (1, 0, 0), 'P': (0, 1, 0),
        'Q': (1, 0, 1), 'R': (0, 1, 1), 'S': (1, 1, 0), 'T': (0, 0, 1),
        'U': (1, 1, 1), 'V': (0, 0, 0), 'W': (1, 0, 0), 'X': (0, 1, 0),
        'Y': (1, 0, 1), 'Z': (0, 1, 1)
    }

    # 1. Calculate the user's profile vector
    user_history = ['A', 'C']
    vec_a = items['A']
    vec_c = items['C']
    profile_vector = tuple( (vec_a[i] + vec_c[i]) / 2 for i in range(3) )

    # 2. Compute cosine similarities for all other items
    recommendation_candidates = {item: vec for item, vec in items.items() if item not in user_history}
    scores = {item: cosine_similarity(profile_vector, vec) for item, vec in recommendation_candidates.items()}

    # 3. Identify best and worst items
    max_score = max(scores.values())
    min_score = min(scores.values())

    best_items = sorted([item for item, score in scores.items() if score == max_score], reverse=True)
    worst_items = [item for item, score in scores.items() if score == min_score]
    
    # Choose F from the worst items based on the prompt's hint
    worst_item = 'F' if 'F' in worst_items else worst_items[0]
    
    recommended_output = ",".join(best_items)
    print(f"Recommended item(s): {recommended_output}")
    print(f"Worst item: {worst_item}")

    # 4. Form an English word
    # The letters are from the best_items (E, M, U) and the worst_item (F)
    word = "FUME"
    print(f"Constructed word: {word}")

    # 5. Encode the word in Base64
    encoded_word = base64.b64encode(word.encode('utf-8')).decode('utf-8')
    print(f"Base64 encoded word: {encoded_word}")

    # 6. Calculate the final vector sum
    word_vectors = [items[char] for char in word]
    final_vector = [0, 0, 0]
    for vec in word_vectors:
        for i in range(3):
            final_vector[i] += vec[i]
    
    # Print the equation with each number
    equation_str = " + ".join(map(str, word_vectors))
    result_str = f"({final_vector[0]},{final_vector[1]},{final_vector[2]})"
    
    print("\nFinal Vector Sum Calculation:")
    print(f"Equation: {equation_str} = {result_str}")
    
if __name__ == '__main__':
    main()
<<<U,M,E F FUME RlVNRQ== (3,3,3)>>>