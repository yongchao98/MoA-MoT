# Jennifer's circle (including Jennifer)
jennifer_males = 6
jennifer_females = 5  # This includes Matilda
total_jennifer_circle = jennifer_males + jennifer_females + 1  # +1 for Jennifer herself

# Thomas's circle (including Thomas)
thomas_males = 4
thomas_females = 1  # This is Matilda
total_thomas_circle = thomas_males + thomas_females + 1  # +1 for Thomas himself

# Since Matilda is in both circles, her female colleagues are:
# From Jennifer's circle: Jennifer herself + all females except Matilda
matilda_females_jennifer = 1 + (jennifer_females - 1)

# From Thomas's circle: No additional females (as she's the only female there)
matilda_females_thomas = 0

# Total female colleagues of Matilda
total_matilda_females = matilda_females_jennifer + matilda_females_thomas

print(f"Matilda has {total_matilda_females} female colleagues")