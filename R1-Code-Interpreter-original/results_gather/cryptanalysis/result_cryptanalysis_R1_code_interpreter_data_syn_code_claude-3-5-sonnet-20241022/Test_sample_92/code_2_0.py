def analyze_number_combinations():
    candidates = [(n, 2) for n in [3, 4, 6]]
    
    def check_against_feedback(n1, n2, guess, feedback):
        g1, g2 = int(guess[0]), int(guess[1])
        
        if feedback == 'both incorrect':
            return n1 != g1 and n2 != g2
        elif feedback == 'both too large':
            return n1 < g1 and n2 < g2
        elif 'one correct position' in feedback:
            correct_count = (n1 == g1) + (n2 == g2)
            return correct_count == 1
        elif 'one wrong position' in feedback:
            wrong_pos_count = (n1 == g2) + (n2 == g1)
            return wrong_pos_count == 1
        elif 'too large' in feedback:
            if g1 > n1 and g2 == n2:
                return True
            if g2 > n2 and g1 == n1:
                return True
        elif 'too small' in feedback:
            if g1 < n1 and g2 == n2:
                return True
            if g2 < n2 and g1 == n1:
                return True
        return True

    feedbacks = [
        ('58TD', 'both incorrect'),
        ('52SK', 'one correct position, one incorrect'),
        ('19UI', 'both incorrect'),
        ('75AC', 'both incorrect'),
        ('86WD', 'one wrong position, one too large'),
        ('72KL', 'one correct position, one too large'),
        ('08LP', 'both incorrect'),
        ('51PB', 'both incorrect'),
        ('96FX', 'one wrong position, one too large'),
        ('02HU', 'one correct position, one too small'),
        ('96FJ', 'one wrong position, one too large'),
        ('12KC', 'one correct position, one too small'),
        ('92TG', 'one correct position, one too large'),
        ('98OJ', 'both too large'),
        ('06ML', 'one wrong position, one too small'),
        ('02VX', 'one correct position, one too small')
    ]

    for n1, n2 in candidates:
        valid = True
        print(f"\nTesting combination: {n1}{n2}")
        for guess, feedback in feedbacks:
            result = check_against_feedback(n1, n2, guess, feedback)
            if not result:
                valid = False
                print(f"Failed at guess {guess}: {feedback}")
                break
        if valid:
            print(f">>> {n1}{n2} satisfies all number conditions <<<")

analyze_number_combinations()